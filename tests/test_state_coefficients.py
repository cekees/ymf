"""Coefficients that are functions of the state: a formula may name the unknowns.

A density ρw (1 - χ) + ρa χ, a conductivity 1 + T², a friction law in h and
q: each is a coefficient whose formula names an unknown. It is substituted
into the equations like any other coefficient, so the ADR form carries its
dependence on the unknowns, and it is not a number, so nothing that wants a
number (a scale, a dimensionless number, boundary data) may lean on it.
"""

from pathlib import Path

import pytest

sympy = pytest.importorskip("sympy")
pytest.importorskip("strictyaml")

from ymf.schema import load_ymf  # noqa: E402
from ymf.symbolic import SymbolicError, adr_problem, state_coefficients  # noqa: E402
from ymf.symbolic.problem import coefficient_numbers  # noqa: E402

EXAMPLES = Path(__file__).parent.parent / "examples"


def _heat(**coefficients):
    """The heat example without its exact solution (it is not exact once κ
    depends on T), with coefficients replaced or added."""
    doc = load_ymf(EXAMPLES / "heat_equation.ymf")
    doc["solution_paths"]["analytical"] = []
    doc["model"]["strong_form"]["coefficients"].update(coefficients)
    return doc


def _equation(problem, component):
    return next(e for e in problem["adr"]["equations"] if e["component"] == component)


def test_a_coefficient_of_the_state_is_the_inlined_expression():
    named = adr_problem(_heat(κ={"formula": "1e-3 (1 + T²)"}))
    inlined = _heat()
    inlined["model"]["strong_form"]["equations"] = ["∂T/∂t = ∇·(1e-3 (1 + T²)∇T)  in Ω"]
    assert named["adr"] == adr_problem(inlined)["adr"]
    diffusion = _equation(named, "T")["diffusion"]["T"]
    assert diffusion["depends_on"] == {"T": "nonlinear"}
    assert diffusion["da"]["T"] == ["0.002*T", "0.002*T"]


def test_a_chain_of_coefficients_reaches_the_state():
    doc = _heat(κ={"formula": "κ0 c"}, c={"formula": "1 + T²"}, κ0={"value": 1.0e-3})
    inlined = _heat()
    inlined["model"]["strong_form"]["equations"] = ["∂T/∂t = ∇·(1e-3 (1 + T²)∇T)  in Ω"]
    assert adr_problem(doc)["adr"] == adr_problem(inlined)["adr"]


def test_a_state_dependent_storage_term_is_a_nonlinear_mass():
    doc = _heat(ρc={"formula": "1 + T"})
    doc["model"]["strong_form"]["equations"] = ["∂(ρc T)/∂t = ∇·(κ∇T)  in Ω"]
    mass = _equation(adr_problem(doc), "T")["mass"]
    assert mass["depends_on"] == {"T": "nonlinear"}


def test_a_coefficient_may_name_a_vector_unknown():
    doc = load_ymf(EXAMPLES / "navier_stokes" / "plane_poiseuille.ymf")
    doc["solution_paths"]["analytical"] = []
    doc["vvuq"] = {}
    strong = doc["model"]["strong_form"]
    strong["coefficients"]["μ"] = {"formula": "1 + v·v"}
    # μΔv is the divergence form only while μ is constant
    strong["equations"][0] = "ρ ∂v/∂t + ∇·(ρ v⊗v) − ∇·(μ∇v) + ∇p = f  in Ω"
    momentum = _equation(adr_problem(doc), "v_0")
    assert set(momentum["diffusion"]["v_0"]["depends_on"]) == {"v_0", "v_1"}


def test_a_cycle_is_still_an_error():
    with pytest.raises(SymbolicError, match="could not be resolved"):
        adr_problem(_heat(κ={"formula": "1e-3 c"}, c={"formula": "κ T"}))


def test_a_coefficient_may_not_take_the_name_of_an_unknown():
    with pytest.raises(SymbolicError, match="name of an unknown"):
        adr_problem(_heat(T={"value": 1.0}))


def test_a_coefficient_may_not_depend_on_a_derivative_yet():
    with pytest.raises(SymbolicError, match="not \\(yet\\) on their derivatives"):
        adr_problem(_heat(κ={"formula": "1e-3 (1 + ∂T/∂x)"}))


@pytest.mark.parametrize("where", ["boundary", "initial", "exact"])
def test_point_data_may_not_lean_on_the_state(where):
    doc = _heat(k={"formula": "300 + T"})
    strong = doc["model"]["strong_form"]
    if where == "boundary":
        strong["boundary_conditions"][0] = {"region": "walls", "variable": "T",
                                            "type": "dirichlet", "formula": "k"}
    elif where == "initial":
        strong["initial_conditions"][0]["formula"] = "k"
    else:
        doc["solution_paths"]["analytical"] = [
            {"name": "wrong", "provenance": "human_specified", "method": "closed_form",
             "solution": {"formula": "T = k"}}]
    with pytest.raises(SymbolicError, match=r"depends on the unknowns \(T\)"):
        adr_problem(doc)


def test_the_exact_solution_of_the_heat_example_is_refused_once_kappa_depends_on_t():
    """Its formula names κ; with κ = κ(T) it is no longer a function of x and t."""
    doc = load_ymf(EXAMPLES / "heat_equation.ymf")
    doc["model"]["strong_form"]["coefficients"]["κ"] = {"formula": "1e-3 (1 + T²)"}
    with pytest.raises(SymbolicError, match="analytical solution"):
        adr_problem(doc)


def test_a_state_coefficient_is_not_a_number():
    spec = {"κ": {"formula": "1e-3 (1 + T²)"}, "κ0": {"value": 2.0}}
    assert coefficient_numbers(spec) == {"κ0": 2.0}


def test_state_coefficients_say_what_they_depend_on():
    strong = _heat(κ={"formula": "κ0 c"}, c={"formula": "1 + T²"}, κ0={"value": 1e-3})
    assert state_coefficients(strong["model"]["strong_form"], 2) == {"κ": ["T"], "c": ["T"]}


def test_a_scale_built_on_a_state_coefficient_is_reported():
    pytest.importorskip("pint")
    from ymf import check_units
    doc = _heat(κ={"formula": "1e-3 (1 + T²)", "units": "m2/s"})
    # the example derives its time scale from κ
    issues = [i for i in check_units(doc).issues
              if i.location == "characteristic_scales[time]"]
    assert len(issues) == 1 and issues[0].severity == "warning"
    assert "coefficients[κ] depends on the unknowns (T)" in issues[0].message
    assert "constant coefficient" in issues[0].message


def test_a_scale_built_on_a_constant_coefficient_is_not_reported():
    pytest.importorskip("pint")
    from ymf import check_units
    doc = _heat(c={"formula": "1 + T²"})
    assert not [i for i in check_units(doc).issues if "depends on the unknowns" in i.message]


def test_a_model_declares_the_coefficient_and_a_case_makes_it_state_dependent(tmp_path):
    (tmp_path / "model.ymf").write_text(
        (EXAMPLES / "heat_equation.ymf").read_text(encoding="utf-8")
        .replace("κ: {value: 1.0e-3, units: m2/s}", "κ: {units: m2/s}")
        .replace('time: {derived: true, formula: "1/(2π²κ)", units: s}',
                 "time: {value: 1.0, units: s}"),
        encoding="utf-8")
    (tmp_path / "case.ymf").write_text(
        "extends: model.ymf\n"
        "remove:\n"
        '  - "solution_paths.analytical[separation_of_variables]"\n'
        "model:\n"
        "  strong_form:\n"
        "    coefficients:\n"
        '      κ: {formula: "1e-3 (1 + T²)"}\n',
        encoding="utf-8")
    assert "κ: {units: m2/s}" in (tmp_path / "model.ymf").read_text(encoding="utf-8")
    doc = load_ymf(tmp_path / "case.ymf")
    assert doc["model"]["strong_form"]["coefficients"]["κ"] == {"units": "m2/s",
                                                                 "formula": "1e-3 (1 + T²)"}
    diffusion = _equation(adr_problem(doc), "T")["diffusion"]["T"]
    assert diffusion["depends_on"] == {"T": "nonlinear"}


# --- the names, kept beside the substituted ADR form --------------------------


def test_a_state_coefficient_keeps_its_name_in_the_plain_data():
    problem = adr_problem(load_ymf(EXAMPLES / "nonlinear_diffusion.ymf"))
    assert problem["coefficients"] == {
        "k": {"code": "u**2 + 1", "depends_on": {"u": "nonlinear"}, "d": {"u": "2*u"},
              "inline": "u**2 + 1"}}


def test_a_problem_without_one_reads_as_before():
    assert "coefficients" not in adr_problem(load_ymf(EXAMPLES / "poisson.ymf"))


def test_every_link_of_a_chain_is_named():
    named = adr_problem(_heat(κ={"formula": "κ0 c"}, c={"formula": "1 + T²"},
                              κ0={"value": 1e-3}))["coefficients"]
    assert set(named) == {"κ", "c"}
    assert named["c"]["d"] == {"T": "2*T"}
    # by name: κ over c, c over T, in that order
    assert named["κ"]["inline"] == "0.001*c" and named["κ"]["uses"] == ["c"]
    assert named["c"]["inline"] == "T**2 + 1"


def test_a_vector_coefficient_of_the_state():
    named = adr_problem(_heat(w={"formula": "(T², T)"}))["coefficients"]
    assert named["w"]["code"] == ["T**2", "T"]
    assert named["w"]["depends_on"] == {"T": "nonlinear"}
    assert named["w"]["d"] == {"T": ["2*T", "1"]}


def test_a_sharp_coefficient_is_named_without_a_derivative():
    named = adr_problem(_heat(χ={"formula": "heaviside(T − 300)"}))["coefficients"]
    assert named["χ"]["sharp"] is True
    assert "d" not in named["χ"]
