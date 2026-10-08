"""The reference graph, and what may be left unsaid.

Phase 3 of the semantics campaign. A YMF document holds competing
approaches to one problem, and the edges between them say which weak form a
discretization implements and which branch a verification plan compares
against. None of those edges was checked: a reference naming nothing still
validated, so a typo in a verification plan was simply ignored.

Two of the edges are single-valued, and where the graph does not branch
they may be omitted -- a discretization in a document with one weak form
can only implement that one. The omission is filled in during composition
and written into the document, so what comes out holds no implicit
bindings.
"""

import pytest

from ymf import load_ymf
from ymf.references import YmfReferenceError

strictyaml = pytest.importorskip("strictyaml")


def spec(weak_forms, discretizations, analytical="", vvuq="", extra=""):
    return """\
kind: model
model:
  name: "p"
  strong_form:
    provenance: llm_derived
    equation_formulation: "−Δu = f"
    unknowns: ["u"]
    domain: "Ω"
    boundary_regions:
      - {name: "all", where: "∂Ω"}
%s%s
solution_paths:
  analytical:%s
  discretizations:%s
%s""" % (("  weak_forms:\n" if weak_forms else ""), weak_forms,
         analytical or " []", discretizations or " []", vvuq)


def weak_form(label):
    return """
    - label: "%s"
      provenance: llm_derived
      derivation: "d"
      solution_spaces: {trial: "H1", test: "H1_0"}
      bilinear: "a"
      linear: "L"
""" % label


def discretization(name, from_weak_form=None):
    line = ('      from_weak_form: "%s"\n' % from_weak_form
            if from_weak_form else "")
    return """
    - name: "%s"
      provenance: llm_derived
%s      solver: {type: linear}""" % (name, line)


def analytical(name):
    return """
    - name: "%s"
      provenance: llm_derived
      method: "closed_form"
      solution: {formula: "0"}""" % name


def write(tmp_path, text, name="p.ymf"):
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


# --------------------------------------------------------------------------
# defaulting where the graph does not branch
# --------------------------------------------------------------------------


def test_one_weak_form_means_a_discretization_need_not_name_it(tmp_path):
    doc = load_ymf(write(tmp_path, spec(weak_form("only"), discretization("d"))))
    assert doc["solution_paths"]["discretizations"][0]["from_weak_form"] == "only"


def test_the_filled_in_reference_is_written_into_the_document(tmp_path):
    """Not resolved on the fly: the document itself carries no implicit binding."""
    doc = load_ymf(write(tmp_path, spec(weak_form("only"), discretization("d"))))
    entry = doc["solution_paths"]["discretizations"][0]
    assert "from_weak_form" in entry


def test_defaulting_is_recorded(tmp_path):
    doc = load_ymf(write(tmp_path, spec(weak_form("only"), discretization("d"))))
    (record,) = doc["composition"]["defaults"]
    assert record["path"] == "solution_paths.discretizations[d].from_weak_form"
    assert record["value"] == "only"
    assert "one weak form" in record["because"]


def test_no_weak_forms_at_all_is_not_an_error(tmp_path):
    """A document may discretize without stating a weak form."""
    doc = load_ymf(write(tmp_path, spec("", discretization("d"))))
    assert "from_weak_form" not in doc["solution_paths"]["discretizations"][0]


def test_two_weak_forms_means_it_has_to_be_said(tmp_path):
    text = spec(weak_form("a") + weak_form("b"), discretization("d"))
    with pytest.raises(YmfReferenceError, match="has 2 weak forms"):
        load_ymf(write(tmp_path, text))


def test_the_error_lists_the_labels_available(tmp_path):
    text = spec(weak_form("a") + weak_form("b"), discretization("d"))
    with pytest.raises(YmfReferenceError, match=r"a, b"):
        load_ymf(write(tmp_path, text))


# --------------------------------------------------------------------------
# "at that point", not "at the end"
# --------------------------------------------------------------------------


def test_a_file_unambiguous_when_written_stays_bound(tmp_path):
    """The reason defaulting happens per file rather than once at the end.

    The model has one weak form, so its discretization is unambiguous. A
    problem extending it adds a second. Resolving against the finished
    document would make the model's omission retroactively ambiguous --
    punishing a file for what a later file did.
    """
    write(tmp_path, spec(weak_form("first"), discretization("d")), "model.ymf")
    variant = write(tmp_path, """\
kind: model
extends: model.ymf
model:
  weak_forms:%s
solution_paths:
  discretizations:%s
""" % (weak_form("second"), discretization("e", "second")), "variant.ymf")
    doc = load_ymf(variant)
    bound = {d["name"]: d["from_weak_form"]
             for d in doc["solution_paths"]["discretizations"]}
    assert bound == {"d": "first", "e": "second"}


def test_removing_the_weak_form_a_binding_points_at_is_caught(tmp_path):
    """Defaulting fills in; resolution checks. remove: can orphan a binding."""
    write(tmp_path, spec(weak_form("only"), discretization("d")), "model.ymf")
    variant = write(tmp_path, """\
kind: model
extends: model.ymf
remove:
  - "model.weak_forms[only]"
""", "variant.ymf")
    with pytest.raises(YmfReferenceError, match="not a weak form label"):
        load_ymf(variant)


# --------------------------------------------------------------------------
# what the graph refuses
# --------------------------------------------------------------------------


def test_a_dangling_weak_form_reference_is_caught(tmp_path):
    text = spec(weak_form("only"), discretization("d", "nonexistent"))
    with pytest.raises(YmfReferenceError, match="not a weak form label"):
        load_ymf(write(tmp_path, text))


def test_an_integer_index_no_longer_names_a_weak_form(tmp_path):
    """Indices were accepted; remove: can renumber them silently."""
    text = spec(weak_form("a") + weak_form("b"), discretization("d", "0"))
    with pytest.raises(YmfReferenceError, match="names '0'"):
        load_ymf(write(tmp_path, text))


def test_an_analytical_entry_may_not_claim_a_weak_form(tmp_path):
    """It solves the strong form; the weak formulation is theory and numerics."""
    text = spec(weak_form("only"), " []",
                analytical("exact").replace(
                    '      method:', '      from_weak_form: "only"\n      method:'))
    with pytest.raises(strictyaml.YAMLValidationError,
                       match="from_weak_form"):
        load_ymf(write(tmp_path, text))


# --------------------------------------------------------------------------
# the verification plan's edges
# --------------------------------------------------------------------------


VVUQ = """
vvuq:
  verification:
    method: "exact_comparison"
%s"""


def test_one_analytical_branch_means_the_reference_need_not_name_it(tmp_path):
    text = spec(weak_form("only"), discretization("d"), analytical("exact"),
                VVUQ % '    reference:\n      type: "analytical"\n')
    doc = load_ymf(write(tmp_path, text))
    assert doc["vvuq"]["verification"]["reference"]["analytical"] == "exact"


def test_the_old_branch_ref_spelling_says_what_it_became(tmp_path):
    text = spec(weak_form("only"), discretization("d"), analytical("exact"),
                VVUQ % '    reference:\n      branch_ref: "exact"\n')
    with pytest.raises(YmfReferenceError, match="renamed to"):
        load_ymf(write(tmp_path, text))


def test_a_reference_to_a_missing_analytical_branch_is_caught(tmp_path):
    text = spec(weak_form("only"), discretization("d"), analytical("exact"),
                VVUQ % '    reference:\n      analytical: "nope"\n')
    with pytest.raises(YmfReferenceError, match="not a analytical solution"):
        load_ymf(write(tmp_path, text))


def test_branches_must_name_discretizations(tmp_path):
    text = spec(weak_form("only"), discretization("d"), " []",
                VVUQ % '    branches: ["d", "ghost"]\n')
    with pytest.raises(YmfReferenceError, match=r"branches\[1\] names 'ghost'"):
        load_ymf(write(tmp_path, text))


def test_expected_rates_keys_must_name_discretizations(tmp_path):
    text = spec(weak_form("only"), discretization("d"), " []",
                VVUQ % '    expected_rates:\n      ghost: {error_L2: 2}\n')
    with pytest.raises(YmfReferenceError, match="expected_rates.ghost"):
        load_ymf(write(tmp_path, text))
