from pathlib import Path

import pytest
import strictyaml

from ymf.schema import load_ymf

EXAMPLES_DIR = Path(__file__).parent.parent / "examples"


def test_load_v01_compatible_example():
    doc = load_ymf(EXAMPLES_DIR / "poisson.ymf")
    assert doc["model"]["name"] == "Poisson unit square"
    # v0.1-style bare string unknowns must still validate
    assert doc["model"]["strong_form"]["unknowns"] == ["u"]


def test_load_v02_enriched_example():
    doc = load_ymf(EXAMPLES_DIR / "poisson_v02_enriched.ymf")
    unknowns = doc["model"]["strong_form"]["unknowns"]
    assert unknowns[0]["name"] == "T"
    assert unknowns[0]["units"] == "K"
    assert unknowns[0]["std_name"] == "land_surface_temperature"
    bc = doc["model"]["strong_form"]["boundary_conditions"][0]
    assert bc["type"] == "dirichlet"
    assert doc["solution_paths"]["discretizations"][0]["bmi_interface"] is True


def test_rejects_unknown_bc_type():
    bad_yaml = """
model:
  name: "bad"
  strong_form:
    provenance: llm_derived
    unknowns: [u]
    equation_formulation: "x"
    strong_form_expression: "x"
    domain: "x"
    boundary_regions:
      - name: "a"
        where: "b"
    boundary_conditions:
      - region: "a"
        variable: "u"
        type: not_a_real_type
solution_paths:
  analytical: []
  discretizations: []
"""
    from ymf.schema import validate_ymf

    with pytest.raises(strictyaml.YAMLValidationError):
        validate_ymf(bad_yaml)


MINIMAL_STRONG_FORM = """
model:
  name: "p"
  strong_form:
{provenance}    unknowns:
      - name: u
{per_unknown}    equation_formulation: "x"
    strong_form_expression: "x"
    domain: "x"
    boundary_regions:
      - name: "a"
        where: "b"
solution_paths:
  analytical: []
  discretizations: []
"""


def _strong_form(provenance="    provenance: llm_derived\n", per_unknown=""):
    return MINIMAL_STRONG_FORM.format(provenance=provenance,
                                      per_unknown=per_unknown)


def test_minimal_strong_form_validates():
    from ymf.schema import validate_ymf

    doc = validate_ymf(_strong_form())
    assert doc["model"]["strong_form"]["provenance"] == "llm_derived"


def test_strong_form_requires_provenance():
    from ymf.schema import validate_ymf

    with pytest.raises(strictyaml.YAMLValidationError):
        validate_ymf(_strong_form(provenance=""))


def test_unknown_provenance_key_is_no_longer_accepted():
    # Renamed to strong_form.provenance; the old key read as "provenance unknown".
    from ymf.schema import validate_ymf

    with pytest.raises(strictyaml.YAMLValidationError):
        validate_ymf(_strong_form(provenance="    unknown_provenance: llm_derived\n"))


def test_an_unknown_may_override_the_strong_form_provenance():
    from ymf.schema import validate_ymf

    doc = validate_ymf(_strong_form(per_unknown="        provenance: human_edited\n"))
    assert doc["model"]["strong_form"]["unknowns"][0]["provenance"] == "human_edited"


def test_an_unknowns_provenance_must_be_a_known_value():
    from ymf.schema import validate_ymf

    with pytest.raises(strictyaml.YAMLValidationError):
        validate_ymf(_strong_form(per_unknown="        provenance: guessed\n"))


def _with_solver(solver_block):
    text = (EXAMPLES_DIR / "poisson.ymf").read_text(encoding="utf-8")
    old = """      solver:
        type: linear
        linear_solver: petsc
        tolerance: 1.0e-12
"""
    assert text.count(old) == 2
    return text.replace(old, solver_block, 1)


def test_petsc_options_are_strings_without_the_dash():
    from ymf.schema import validate_ymf

    doc = validate_ymf(_with_solver("""      solver:
        type: linear
        linear_solver: petsc
        petsc_options: {ksp_type: cg, pc_type: gamg, ksp_rtol: "1e-10", ksp_view: ""}
"""))
    options = doc["solution_paths"]["discretizations"][0]["solver"]["petsc_options"]
    assert options == {"ksp_type": "cg", "pc_type": "gamg",
                       "ksp_rtol": "1e-10", "ksp_view": ""}


def test_petsc_options_must_be_a_mapping():
    from ymf.schema import validate_ymf

    with pytest.raises(strictyaml.YAMLValidationError):
        validate_ymf(_with_solver("""      solver:
        type: linear
        linear_solver: petsc
        petsc_options: [ksp_type, cg]
"""))
