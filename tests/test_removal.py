"""Composition may remove, not only add and change.

Phase 2 of the semantics campaign. A problem could state more than the
model it extends, and could change what the model said, but could not drop
anything: taking the ``incompressible`` assumption out of a Navier-Stokes
model, to solve it as compressible, meant writing a new model file.

The motivating moves are real ones. Incompressible flow solved as
compressible, or inviscid Euler as viscous; and in the other direction an
unknown field dropped, or replaced by an analytical expression.

Every removal is recorded in ``composition.removals``, beside the
``overrides`` that record every change, so a composed document still says
exactly how it came to be.
"""

import pytest

from ymf import load_ymf
from ymf.compose import YmfRemovalError

strictyaml = pytest.importorskip("strictyaml")


MODEL = """\
ymf_state: abstract
model:
  name: "flow"
  physical_model:
    provenance: llm_derived
    processes: [momentum_transport]
    assumptions: [incompressible, viscous, steady]
  strong_form:
    provenance: llm_derived
    equation_formulation: "ρ(v·∇)v − μΔv + ∇p = f"
    unknowns: ["v", "p"]
    domain: "Ω"
    boundary_regions:
      - {name: "walls", where: "∂Ω"}
    coefficients:
      ρ: {value: 1.0, units: "kg/m3"}
      μ: {value: 1.0, units: "Pa*s"}
  weak_forms:
    - label: "galerkin"
      provenance: llm_derived
      derivation: "multiply by w and integrate by parts"
      solution_spaces: {trial: "H1", test: "H1_0"}
      bilinear: "a(v,w)"
      linear: "L(w)"
    - label: "stabilized"
      provenance: llm_derived
      derivation: "add SUPG terms"
      solution_spaces: {trial: "H1", test: "H1_0"}
      bilinear: "a(v,w) + s(v,w)"
      linear: "L(w)"
"""


def write(tmp_path, variant_body):
    (tmp_path / "model.ymf").write_text(MODEL, encoding="utf-8")
    path = tmp_path / "variant.ymf"
    path.write_text("ymf_state: abstract\nextends: model.ymf\n" + variant_body,
                    encoding="utf-8")
    return path


# --------------------------------------------------------------------------
# the motivating moves
# --------------------------------------------------------------------------


def test_an_assumption_can_be_dropped(tmp_path):
    """Incompressible flow, solved as compressible."""
    doc = load_ymf(write(tmp_path, """\
remove:
  - "model.physical_model.assumptions[incompressible]"
"""))
    assert doc["model"]["physical_model"]["assumptions"] == ["viscous", "steady"]


def test_dropping_and_adding_in_one_step(tmp_path):
    """The merge unions the lists first, then the removal takes the one named."""
    doc = load_ymf(write(tmp_path, """\
remove:
  - "model.physical_model.assumptions[incompressible]"
model:
  physical_model:
    assumptions: [compressible]
"""))
    assert doc["model"]["physical_model"]["assumptions"] == [
        "viscous", "steady", "compressible"]


def test_an_unknown_field_can_be_dropped(tmp_path):
    doc = load_ymf(write(tmp_path, """\
remove:
  - "model.strong_form.unknowns[p]"
"""))
    assert doc["model"]["strong_form"]["unknowns"] == ["v"]


def test_a_coefficient_can_be_dropped(tmp_path):
    """Euler: the same model without its viscosity."""
    doc = load_ymf(write(tmp_path, """\
remove:
  - "model.strong_form.coefficients.μ"
"""))
    assert set(doc["model"]["strong_form"]["coefficients"]) == {"ρ"}


def test_a_labelled_list_entry_can_be_dropped(tmp_path):
    doc = load_ymf(write(tmp_path, """\
remove:
  - "model.weak_forms[stabilized]"
"""))
    assert [w["label"] for w in doc["model"]["weak_forms"]] == ["galerkin"]


def test_several_paths_in_one_file(tmp_path):
    doc = load_ymf(write(tmp_path, """\
remove:
  - "model.physical_model.assumptions[viscous]"
  - "model.strong_form.coefficients.μ"
"""))
    assert doc["model"]["physical_model"]["assumptions"] == [
        "incompressible", "steady"]
    assert set(doc["model"]["strong_form"]["coefficients"]) == {"ρ"}


# --------------------------------------------------------------------------
# the record
# --------------------------------------------------------------------------


def test_a_removal_is_recorded_with_what_it_took_and_who_took_it(tmp_path):
    doc = load_ymf(write(tmp_path, """\
remove:
  - "model.strong_form.coefficients.μ"
"""))
    assert doc["composition"]["removals"] == [{
        "path": "model.strong_form.coefficients.μ",
        "was": {"value": 1.0, "units": "Pa*s"},
        "removed_by": "variant.ymf",
    }]


def test_removals_sit_beside_overrides_not_inside_them(tmp_path):
    """A consumer reading overrides expects every entry to have a `now`."""
    doc = load_ymf(write(tmp_path, """\
remove:
  - "model.strong_form.coefficients.μ"
model:
  strong_form:
    coefficients:
      ρ: {value: 2.0}
"""))
    composition = doc["composition"]
    assert all("now" in o for o in composition["overrides"])
    assert [r["path"] for r in composition["removals"]] == [
        "model.strong_form.coefficients.μ"]


def test_a_document_that_removes_nothing_records_nothing(tmp_path):
    doc = load_ymf(write(tmp_path, """\
model:
  strong_form:
    coefficients:
      ρ: {value: 2.0}
"""))
    assert "removals" not in doc["composition"]


# --------------------------------------------------------------------------
# a path that names nothing
# --------------------------------------------------------------------------


@pytest.mark.parametrize("path,fragment", [
    ("model.strong_form.coefficients.κ", "nothing to remove"),
    ("model.physical_model.assumptions[turbulent]", "no entry 'turbulent'"),
    ("model.weak_forms[petrov]", "no entry 'petrov'"),
    ("model.nonexistent.deeper", "nothing at"),
    ("model.strong_form.domain[x]", "is not a list"),
])
def test_removing_what_is_not_there_is_an_error(tmp_path, path, fragment):
    """A stale path means the document moved under the file removing from it.

    Composing the unremoved thing silently is how a variant ends up stating
    the opposite of what its author wrote.
    """
    with pytest.raises(YmfRemovalError, match=fragment):
        load_ymf(write(tmp_path, 'remove:\n  - "%s"\n' % path))


def test_a_malformed_path_is_refused(tmp_path):
    with pytest.raises(YmfRemovalError, match="unclosed"):
        load_ymf(write(tmp_path, 'remove:\n  - "model.weak_forms[galerkin"\n'))
