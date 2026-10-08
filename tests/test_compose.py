"""Composing a YMF document from several files with ``extends``."""

from pathlib import Path

import pytest
import strictyaml

from ymf.compose import YmfCompositionError
from ymf.schema import load_ymf

NAVIER_STOKES = Path(__file__).parent.parent / "examples" / "navier_stokes"

MODEL = """
kind: model
Problem:
  physical_model:
    provenance: human_specified
    assumptions: [newtonian]
  strong_form:
    provenance: human_specified
    unknowns:
      - name: u
        units: m/s
    equation_formulation: "model"
    strong_form_expression: "x"
    coefficients:
      mu: {units: Pa*s}
"""

PROBLEM = """
extends: model.ymf
Problem:
  name: "problem"
  physical_model:
    provenance: human_specified
    assumptions: [steady]
  strong_form:
    domain: "Ω"
    boundary_regions:
      - name: "b"
        geometry: "∂Ω"
    coefficients:
      mu: {value: 1.0}
solution_paths:
  analytical: []
  discretizations: []
"""


def write(tmp_path, **files):
    for name, text in files.items():
        (tmp_path / (name + ".ymf")).write_text(text, encoding="utf-8")
    return tmp_path


# --- the examples ----------------------------------------------------------


def test_the_model_loads_as_a_partial_document():
    doc = load_ymf(NAVIER_STOKES / "navier_stokes_model.ymf")
    assert doc["kind"] == "model"
    assert "name" not in doc["Problem"]
    assert "solution_paths" not in doc
    assert "composition" not in doc


def test_couette_composes_the_model_with_its_own_additions():
    doc = load_ymf(NAVIER_STOKES / "planar_couette.ymf")
    strong = doc["Problem"]["strong_form"]
    # from the model
    assert strong["equations"][0].startswith("ρ ∂v/∂t")
    assert [u["name"] for u in strong["unknowns"]] == ["v", "p"]
    # merged: the model's units, the problem's value
    # a number, not the string '0.1': coefficients are typed in the schema
    assert strong["coefficients"]["μ"] == {"units": "Pa*s", "value": 0.1}
    # string lists are unioned in order
    assert doc["Problem"]["physical_model"]["assumptions"] == [
        "newtonian", "constant_density", "isothermal", "steady", "fully_developed"]
    assert {k: v for k, v in doc["composition"].items() if k != "defaults"} == {
        "sources": ["navier_stokes_model.ymf", "planar_couette.ymf"]}


def test_a_parameter_variant_records_exactly_what_it_overrides():
    doc = load_ymf(NAVIER_STOKES / "plane_poiseuille_re100.ymf")
    assert doc["composition"]["sources"] == [
        "navier_stokes_model.ymf", "plane_poiseuille.ymf", "plane_poiseuille_re100.ymf"]
    by_path = {o["path"]: o for o in doc["composition"]["overrides"]}
    # The parent derives its scales from μ and G, so those two and the
    # name are all the variant has to change.
    assert set(by_path) == {
        "Problem.name",
        "Problem.strong_form.coefficients.μ.value",
        "Problem.strong_form.coefficients.G.value",
    }
    assert by_path["Problem.strong_form.coefficients.μ.value"] == {
        "path": "Problem.strong_form.coefficients.μ.value",
        "was": 1.0, "now": 0.01, "set_by": "plane_poiseuille_re100.ymf"}


def test_the_composed_variant_resolves_its_reynolds_number():
    pytest.importorskip("pint")   # resolving it is ymf.units' job
    from ymf.units import compute_dimensionless_numbers

    for name, re in [("planar_couette", 10.0), ("plane_poiseuille", 1.0),
                     ("plane_poiseuille_re100", 100.0)]:
        doc = load_ymf(NAVIER_STOKES / (name + ".ymf"))
        assert compute_dimensionless_numbers(doc["Problem"])["Re"] == pytest.approx(re)


# --- merge rules -----------------------------------------------------------


def test_list_entries_merge_by_name_and_new_ones_append(tmp_path):
    write(tmp_path, model=MODEL, problem=PROBLEM.replace(
        "    coefficients:", """    unknowns:
      - name: u
        std_name: x__velocity
      - name: p
        units: Pa
    coefficients:"""))
    doc = load_ymf(tmp_path / "problem.ymf")
    assert doc["Problem"]["strong_form"]["unknowns"] == [
        {"name": "u", "units": "m/s", "std_name": "x__velocity"},
        {"name": "p", "units": "Pa"},
    ]
    assert "overrides" not in doc["composition"]


def test_changing_a_list_entry_field_is_an_override_named_by_identity(tmp_path):
    write(tmp_path, model=MODEL, problem=PROBLEM.replace(
        "    coefficients:", """    unknowns:
      - name: u
        units: km/h
    coefficients:"""))
    doc = load_ymf(tmp_path / "problem.ymf")
    assert doc["composition"]["overrides"] == [{
        "path": "Problem.strong_form.unknowns[u].units",
        "was": "m/s", "now": "km/h", "set_by": "problem.ymf"}]


def test_file_level_keys_are_not_merged(tmp_path):
    write(tmp_path, model=MODEL, problem=PROBLEM)
    doc = load_ymf(tmp_path / "problem.ymf")
    assert "extends" not in doc and "kind" not in doc


def test_a_model_may_extend_a_model(tmp_path):
    write(tmp_path, model=MODEL, refined="""
kind: model
extends: model.ymf
Problem:
  physical_model:
    provenance: human_specified
    assumptions: [isothermal]
""")
    doc = load_ymf(tmp_path / "refined.ymf")
    assert doc["Problem"]["physical_model"]["assumptions"] == ["newtonian", "isothermal"]
    assert doc["kind"] == "model"
    assert doc["composition"]["sources"] == ["model.ymf", "refined.ymf"]


# --- validation and errors -------------------------------------------------


def test_an_incomplete_composition_is_rejected_naming_the_chain(tmp_path):
    write(tmp_path, model=MODEL, problem=PROBLEM.replace('  name: "problem"\n', ""))
    with pytest.raises(strictyaml.YAMLValidationError, match=r"composed from model\.ymf, problem\.ymf"):
        load_ymf(tmp_path / "problem.ymf")


def test_a_model_on_its_own_is_not_a_problem(tmp_path):
    write(tmp_path, model=MODEL.replace("kind: model\n", ""))
    with pytest.raises(strictyaml.YAMLValidationError):
        load_ymf(tmp_path / "model.ymf")


def test_an_error_in_a_parent_is_reported_against_the_parent(tmp_path):
    write(tmp_path, model=MODEL.replace("provenance: human_specified\n    assumptions",
                                        "provenance: guessed\n    assumptions"),
          problem=PROBLEM)
    with pytest.raises(strictyaml.YAMLValidationError, match=r"model\.ymf"):
        load_ymf(tmp_path / "problem.ymf")


def test_an_error_in_a_single_file_is_reported_against_it(tmp_path):
    path = tmp_path / "bad.ymf"
    path.write_text("Problem:\n  name: x\n  nonsense: 1\nsolution_paths:\n"
                    "  analytical: []\n  discretizations: []\n", encoding="utf-8")
    with pytest.raises(strictyaml.YAMLValidationError, match=r"bad\.ymf"):
        load_ymf(path)


def test_a_missing_parent_is_named(tmp_path):
    write(tmp_path, problem=PROBLEM)
    with pytest.raises(YmfCompositionError, match=r"extends .*model\.ymf, which does not exist"):
        load_ymf(tmp_path / "problem.ymf")


def test_an_extends_cycle_is_rejected(tmp_path):
    write(tmp_path, a="kind: model\nextends: b.ymf\n", b="kind: model\nextends: a.ymf\n")
    with pytest.raises(YmfCompositionError, match="cycle"):
        load_ymf(tmp_path / "a.ymf")


def test_a_typed_override_keeps_its_type(tmp_path):
    write(tmp_path, model=MODEL, problem=PROBLEM.replace(
        "  strong_form:", "  characteristic_scales:\n    length: {value: 1.0}\n  strong_form:"),
        variant="""
extends: problem.ymf
Problem:
  characteristic_scales:
    length: {value: 2.5}
""")
    doc = load_ymf(tmp_path / "variant.ymf")
    assert doc["composition"]["overrides"] == [{
        "path": "Problem.characteristic_scales.length.value",
        "was": 1.0, "now": 2.5, "set_by": "variant.ymf"}]


# ---------------------------------------------------------------------------
# attributing a composed document's failure to a file
# ---------------------------------------------------------------------------
#
# A composed problem has no file of its own. strictyaml could only point at a
# line of a document that exists nowhere on disk, labelled with the chain, so
# the reader was left to work out which of several files to open. The failure
# that only a composition can have is a required key that no file in the
# chain sets -- each file is valid on its own against the partial schema, and
# only the whole is checked in full.


BASE = """\
kind: model
Problem:
  name: "flow"
  strong_form:
    provenance: llm_derived
    equation_formulation: "f"
    unknowns: ["v"]
    domain: "Ω"
    boundary_regions:
      - {name: "w", geometry: "∂Ω"}
solution_paths:
  analytical: []
  discretizations: []
"""


def test_a_missing_required_key_names_the_files_that_write_its_mapping(tmp_path):
    write(tmp_path,
          model=BASE.replace('    equation_formulation: "f"\n', ""),
          problem="""
extends: model.ymf
Problem:
  strong_form:
    coefficients:
      ρ: {value: 1.0}
""")
    with pytest.raises(strictyaml.YAMLValidationError) as caught:
        load_ymf(tmp_path / "problem.ymf")
    note = str(caught.value)
    assert "no file in the chain sets equation_formulation under Problem.strong_form" in note
    # both files write that mapping, named in chain order
    assert "Problem.strong_form is written by model.ymf, problem.ymf" in note


def test_an_incomplete_new_entry_names_the_file_that_added_it(tmp_path):
    """The useful case: one file added it, so there is one file to open."""
    write(tmp_path, model=BASE, problem="""
extends: model.ymf
solution_paths:
  discretizations:
    - name: "new_path"
      provenance: llm_derived
""")
    with pytest.raises(strictyaml.YAMLValidationError) as caught:
        load_ymf(tmp_path / "problem.ymf")
    note = str(caught.value)
    assert "solution_paths.discretizations[new_path] is written by problem.ymf" in note
    assert "model.ymf" not in note.split("is written by")[-1]


def test_the_attributed_path_is_the_one_remove_and_overrides_use(tmp_path):
    """So a path out of an error can be pasted into remove:, and vice versa."""
    write(tmp_path, model=BASE, problem="""
extends: model.ymf
solution_paths:
  discretizations:
    - name: "new_path"
      provenance: llm_derived
""")
    with pytest.raises(strictyaml.YAMLValidationError) as caught:
        load_ymf(tmp_path / "problem.ymf")
    # "[new_path]", not ".name" -- the entry, not the key carrying its identity
    assert "discretizations[new_path]" in str(caught.value)


def test_a_file_that_is_wrong_on_its_own_still_cites_itself(tmp_path):
    """Unchanged: per-file errors never needed attributing.

    Each file is validated as written before anything is merged, so an error
    a single file holds is reported against that file's own line. Only what
    no file holds alone reaches the composed check.
    """
    write(tmp_path, model=BASE, problem="""
extends: model.ymf
Problem:
  strong_form:
    provenance: not_a_provenance
""")
    with pytest.raises(strictyaml.YAMLValidationError) as caught:
        load_ymf(tmp_path / "problem.ymf")
    message = str(caught.value)
    assert "problem.ymf" in message
    assert "composed from" not in message
