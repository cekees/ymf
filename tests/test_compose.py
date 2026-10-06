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
extends: model.yaml
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
        (tmp_path / (name + ".yaml")).write_text(text, encoding="utf-8")
    return tmp_path


# --- the examples ----------------------------------------------------------


def test_the_model_loads_as_a_partial_document():
    doc = load_ymf(NAVIER_STOKES / "navier_stokes_model.yaml")
    assert doc["kind"] == "model"
    assert "name" not in doc["Problem"]
    assert "solution_paths" not in doc
    assert "composition" not in doc


def test_couette_composes_the_model_with_its_own_additions():
    doc = load_ymf(NAVIER_STOKES / "planar_couette.yaml")
    strong = doc["Problem"]["strong_form"]
    # from the model
    assert strong["strong_form_expression"].startswith("ρ (∂v/∂t")
    assert [u["name"] for u in strong["unknowns"]] == ["v", "p"]
    # merged: the model's units, the problem's value
    assert strong["coefficients"]["mu"] == {"units": "Pa*s", "value": "0.1"}
    # string lists are unioned in order
    assert doc["Problem"]["physical_model"]["assumptions"] == [
        "newtonian", "constant_density", "isothermal", "steady", "fully_developed"]
    assert doc["composition"] == {
        "sources": ["navier_stokes_model.yaml", "planar_couette.yaml"]}


def test_a_parameter_variant_records_exactly_what_it_overrides():
    doc = load_ymf(NAVIER_STOKES / "plane_poiseuille_re100.yaml")
    assert doc["composition"]["sources"] == [
        "navier_stokes_model.yaml", "plane_poiseuille.yaml", "plane_poiseuille_re100.yaml"]
    by_path = {o["path"]: o for o in doc["composition"]["overrides"]}
    # The parent derives its scales from mu and G, so those two and the
    # name are all the variant has to change.
    assert set(by_path) == {
        "Problem.name",
        "Problem.strong_form.coefficients.mu.value",
        "Problem.strong_form.coefficients.G.value",
    }
    assert by_path["Problem.strong_form.coefficients.mu.value"] == {
        "path": "Problem.strong_form.coefficients.mu.value",
        "was": "1.0", "now": "0.01", "set_by": "plane_poiseuille_re100.yaml"}


def test_the_composed_variant_resolves_its_reynolds_number():
    from ymf.units import compute_dimensionless_numbers

    for name, re in [("planar_couette", 10.0), ("plane_poiseuille", 1.0),
                     ("plane_poiseuille_re100", 100.0)]:
        doc = load_ymf(NAVIER_STOKES / (name + ".yaml"))
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
    doc = load_ymf(tmp_path / "problem.yaml")
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
    doc = load_ymf(tmp_path / "problem.yaml")
    assert doc["composition"]["overrides"] == [{
        "path": "Problem.strong_form.unknowns[u].units",
        "was": "m/s", "now": "km/h", "set_by": "problem.yaml"}]


def test_file_level_keys_are_not_merged(tmp_path):
    write(tmp_path, model=MODEL, problem=PROBLEM)
    doc = load_ymf(tmp_path / "problem.yaml")
    assert "extends" not in doc and "kind" not in doc


def test_a_model_may_extend_a_model(tmp_path):
    write(tmp_path, model=MODEL, refined="""
kind: model
extends: model.yaml
Problem:
  physical_model:
    provenance: human_specified
    assumptions: [isothermal]
""")
    doc = load_ymf(tmp_path / "refined.yaml")
    assert doc["Problem"]["physical_model"]["assumptions"] == ["newtonian", "isothermal"]
    assert doc["kind"] == "model"
    assert doc["composition"]["sources"] == ["model.yaml", "refined.yaml"]


# --- validation and errors -------------------------------------------------


def test_an_incomplete_composition_is_rejected_naming_the_chain(tmp_path):
    write(tmp_path, model=MODEL, problem=PROBLEM.replace('  name: "problem"\n', ""))
    with pytest.raises(strictyaml.YAMLValidationError, match=r"composed from model\.yaml, problem\.yaml"):
        load_ymf(tmp_path / "problem.yaml")


def test_a_model_on_its_own_is_not_a_problem(tmp_path):
    write(tmp_path, model=MODEL.replace("kind: model\n", ""))
    with pytest.raises(strictyaml.YAMLValidationError):
        load_ymf(tmp_path / "model.yaml")


def test_an_error_in_a_parent_is_reported_against_the_parent(tmp_path):
    write(tmp_path, model=MODEL.replace("provenance: human_specified\n    assumptions",
                                        "provenance: guessed\n    assumptions"),
          problem=PROBLEM)
    with pytest.raises(strictyaml.YAMLValidationError, match=r"model\.yaml"):
        load_ymf(tmp_path / "problem.yaml")


def test_an_error_in_a_single_file_is_reported_against_it(tmp_path):
    path = tmp_path / "bad.yaml"
    path.write_text("Problem:\n  name: x\n  nonsense: 1\nsolution_paths:\n"
                    "  analytical: []\n  discretizations: []\n", encoding="utf-8")
    with pytest.raises(strictyaml.YAMLValidationError, match=r"bad\.yaml"):
        load_ymf(path)


def test_a_missing_parent_is_named(tmp_path):
    write(tmp_path, problem=PROBLEM)
    with pytest.raises(YmfCompositionError, match=r"extends .*model\.yaml, which does not exist"):
        load_ymf(tmp_path / "problem.yaml")


def test_an_extends_cycle_is_rejected(tmp_path):
    write(tmp_path, a="kind: model\nextends: b.yaml\n", b="kind: model\nextends: a.yaml\n")
    with pytest.raises(YmfCompositionError, match="cycle"):
        load_ymf(tmp_path / "a.yaml")


def test_a_typed_override_keeps_its_type(tmp_path):
    write(tmp_path, model=MODEL, problem=PROBLEM.replace(
        "  strong_form:", "  characteristic_scales:\n    length: {value: 1.0}\n  strong_form:"),
        variant="""
extends: problem.yaml
Problem:
  characteristic_scales:
    length: {value: 2.5}
""")
    doc = load_ymf(tmp_path / "variant.yaml")
    assert doc["composition"]["overrides"] == [{
        "path": "Problem.characteristic_scales.length.value",
        "was": 1.0, "now": 2.5, "set_by": "variant.yaml"}]
