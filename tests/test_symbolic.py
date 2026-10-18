"""ymf.symbolic: strong form -> advection-diffusion-reaction form."""

from pathlib import Path

import pytest
import yaml

sympy = pytest.importorskip("sympy")
pytest.importorskip("strictyaml")

from ymf.schema import load_ymf  # noqa: E402
from ymf.symbolic import (Space, SymbolicError, adr_form, adr_problem,  # noqa: E402
                          parse_equation, parse_expression)

EXAMPLES = Path(__file__).parent.parent / "examples"


def classify(text, unknowns, coefficients=None, dim=2):
    space = Space(dim, unknowns, coefficients)
    residuals = []
    for line in text if isinstance(text, list) else [text]:
        residuals.extend(parse_equation(line, space))
    return adr_form(residuals, space)


# --- classification ---------------------------------------------------------


def test_poisson_is_diagonal_diffusion_and_a_constant_reaction():
    out = classify("-div(grad(u)) = 2*pi**2*sin(pi*x)*sin(pi*y)", {"u": 0})
    (eq,) = out["equations"]
    assert set(eq) == {"component", "diffusion", "reaction"}
    d = eq["diffusion"]["u"]
    assert (d["rowptr"], d["colind"], d["a"]) == ([0, 1, 2], [0, 1], ["1", "1"])
    assert d["depends_on"] == {}            # constant in the unknowns
    assert eq["reaction"]["r"] == "-2*numpy.pi**2*numpy.sin(numpy.pi*x)*numpy.sin(numpy.pi*y)"
    assert eq["reaction"]["depends_on"] == {}


def test_a_variable_coefficient_is_differentiated_not_pulled_out():
    # div(x*grad(u)) = x*lap(u) + du/dx: the coefficient stays inside the
    # flux, so the form is diffusion with a = x, not x * (a second derivative)
    out = classify("-div(x*grad(u)) = 0", {"u": 0})
    assert out["equations"][0]["diffusion"]["u"]["a"] == ["x", "x"]


def test_transient_nonlinear_advection_diffusion():
    out = classify("dt(u) + div(u**2*(1, 0)) - div(grad(u)) = 0", {"u": 0})
    eq = out["equations"][0]
    assert eq["mass"] == {"m": "u", "dm": {"u": "1"}, "depends_on": {"u": "linear"}}
    assert eq["advection"]["f"] == ["u**2", "0"]
    assert eq["advection"]["df"] == {"u": ["2*u", "0"]}
    assert eq["advection"]["depends_on"] == {"u": "nonlinear"}


def test_navier_stokes_matches_proteus_structure():
    out = classify(["rho*dt(v) + div(rho*outer(v, v)) - div(mu*grad(v)) + grad(p) = f",
                    "div(v) = 0"],
                   {"v": 1, "p": 0},
                   {"rho": sympy.Float(1), "mu": sympy.Float(0.1), "f": sympy.Matrix([0, 0])})
    assert out["components"] == ["v_0", "v_1", "p"]
    v0, v1, p = out["equations"]
    assert v0["mass"]["depends_on"] == {"v_0": "linear"}
    assert v0["advection"]["depends_on"] == {"v_0": "nonlinear", "v_1": "nonlinear", "p": "linear"}
    assert v0["diffusion"]["v_0"]["a"] == ["0.1", "0.1"]
    # continuity: the velocity is the advective flux of the pressure equation
    assert set(p) == {"component", "advection"}
    assert p["advection"]["f"] == ["v_0", "v_1"]
    assert p["advection"]["depends_on"] == {"v_0": "linear", "v_1": "linear"}


def test_a_gradient_hamiltonian_and_its_derivatives():
    # an eikonal-like term, not in divergence form
    out = classify("dt(u) + dx(u)**2 + dy(u)**2 = 1", {"u": 0})
    h = out["equations"][0]["hamiltonian"]
    assert h["H"] == "grad_u_0**2 + grad_u_1**2"
    assert h["dH"] == {"u": ["2*grad_u_0", "2*grad_u_1"]}
    assert h["depends_on"] == {"u": "nonlinear"}


def test_a_pressure_gradient_can_be_a_hamiltonian_instead_of_a_flux():
    doc = load_ymf(EXAMPLES / "kovasznay_flow.ymf")
    flux = {e["component"]: e for e in adr_problem(doc)["adr"]["equations"]}
    ham = {e["component"]: e for e in
           adr_problem(doc, hamiltonian_gradients=("p",))["adr"]["equations"]}
    assert flux["v_0"]["advection"]["f"][0] == "p + 1.0*v_0**2" and "hamiltonian" not in flux["v_0"]
    assert ham["v_0"]["advection"]["f"] == ["1.0*v_0**2", "1.0*v_0*v_1"]
    assert ham["v_0"]["hamiltonian"] == {"H": "grad_p_0", "dH": {"p": ["1", "0"]},
                                         "depends_on": {"p": "linear"}}
    assert ham["p"] == flux["p"]          # continuity is untouched


def test_an_unknown_seen_only_through_its_gradient_and_never_fixed_is_up_to_a_constant():
    assert adr_problem(load_ymf(EXAMPLES / "kovasznay_flow.ymf"))["up_to_constant"] == ["p"]
    # a pressure datum fixes it; Poisson's u has Dirichlet data
    assert "up_to_constant" not in adr_problem(load_ymf(EXAMPLES / "navier_stokes" / "plane_poiseuille.ymf"))
    assert "up_to_constant" not in adr_problem(load_ymf(EXAMPLES / "poisson.ymf"))


def test_the_output_is_plain_data():
    out = classify(["rho*dt(v) + div(rho*outer(v, v)) - div(mu*grad(v)) + grad(p) = f",
                    "div(v) = 0"],
                   {"v": 1, "p": 0},
                   {"rho": sympy.Float(1), "mu": sympy.Float(1), "f": sympy.Matrix([1, 0])})
    assert yaml.safe_load(yaml.safe_dump(out)) == out


# --- what it refuses, and how it says so -------------------------------------


@pytest.mark.parametrize("equations, unknowns, message", [
    ("-div(grad(u)) = g", {"u": 0}, r"unknown name 'g'"),
    # (dx(dx(u)) alone *is* a divergence: dx(1*dx(u)), diffusion a = diag(-1, 0))
    ("u*dx(dx(u)) = 0", {"u": 0}, r"second derivatives .* not in divergence form"),
    # the convective form: the Hamiltonian would depend on v, not only grad v
    (["dot(v, grad(v)) - div(grad(v)) + grad(p) = 0", "div(v) = 0"], {"v": 1, "p": 0},
     r"depends on the unknowns themselves"),
    (["-div(grad(u)) = 0"], {"u": 0, "w": 0}, r"1 scalar equations for 2 scalar unknowns"),
    ("u*dt(u) = 0", {"u": 0}, r"time derivative .* not of the form c\*dt\(e\)"),
])
def test_unsupported_forms_are_refused_with_a_reason(equations, unknowns, message):
    with pytest.raises(SymbolicError, match=message):
        classify(equations, unknowns)


def test_parse_errors_name_the_text():
    with pytest.raises(SymbolicError, match=r"expected a number, a name or '\(' at position 5 in 'div\(\('"):
        parse_expression("div((", Space(2, {"u": 0}))


# --- whole specs ---------------------------------------------------------------


def test_the_poisson_example_becomes_a_runnable_problem():
    problem = adr_problem(load_ymf(EXAMPLES / "poisson.ymf"))
    assert problem["domain"] == {"lower": [0.0, 0.0], "upper": [1.0, 1.0]}
    assert problem["dirichlet"] == [{"component": "u", "region": "all_boundaries",
                                     "where": None, "value": "0.0"}]
    assert problem["exact"] == {"u": "(1/2)*numpy.sin(numpy.pi*x)*numpy.sin(numpy.pi*y)"}
    assert problem["time"] is None


def _poisson_with(bc):
    """The Poisson example with its one boundary condition replaced."""
    doc = load_ymf(EXAMPLES / "poisson.ymf")
    doc["model"]["strong_form"]["boundary_conditions"] = [bc]
    return doc


def test_a_neumann_condition_becomes_an_outward_flux():
    """Its value is q.n with q = -a grad u, which is the sign Proteus wants."""
    problem = adr_problem(_poisson_with(
        {"region": "all_boundaries", "variable": "u", "type": "neumann",
         "formula": "-2*x"}))
    assert problem["dirichlet"] == []
    assert problem["neumann"] == [{"component": "u", "region": "all_boundaries",
                                   "where": None, "value": "-2*x"}]


def test_a_neumann_only_poisson_is_determined_up_to_a_constant():
    """No Dirichlet condition fixes u anywhere, so the constant is free."""
    problem = adr_problem(_poisson_with(
        {"region": "all_boundaries", "variable": "u", "type": "neumann",
         "value": 0.0}))
    assert problem.get("up_to_constant") == ["u"]


def test_a_discharge_hydrograph_is_a_tabulated_neumann():
    problem = adr_problem(_poisson_with(
        {"region": "all_boundaries", "variable": "u", "type": "neumann",
         "series": [[0.0, -2.0], [1.0, -4.0]]}))
    assert problem["neumann"] == [{"component": "u", "region": "all_boundaries",
                                   "where": None,
                                   "series": [[0.0, -2.0], [1.0, -4.0]]}]


def test_a_stage_hydrograph_is_a_tabulated_dirichlet():
    """The case naming "hydrograph" as a type made impossible to write.

    A gauge reading against time is a value, not a flux, so it belongs under
    dirichlet. Supplying a value by table is orthogonal to the kind of
    condition it is.
    """
    problem = adr_problem(_poisson_with(
        {"region": "all_boundaries", "variable": "u", "type": "dirichlet",
         "series": [[0.0, 1.0], [1.0, 2.0]]}))
    assert problem["dirichlet"] == [{"component": "u", "region": "all_boundaries",
                                     "where": None,
                                     "series": [[0.0, 1.0], [1.0, 2.0]]}]
    assert problem["neumann"] == []


def test_hydrograph_is_not_a_condition_type():
    with pytest.raises(Exception, match="hydrograph"):
        adr_problem(_poisson_with(
            {"region": "all_boundaries", "variable": "u", "type": "hydrograph",
             "series": [[0.0, -2.0], [1.0, -4.0]]}))


@pytest.mark.parametrize("series,message", [
    ([[0.0, -2.0]], "one sample is a constant"),
    ([[1.0, -2.0], [0.0, -4.0]], "times must increase"),
    ([[0.0, -2.0], [0.0, -4.0]], "times must increase"),
    ([[0.0, -2.0, 7.0], [1.0, -4.0, 9.0]], "each sample is"),
])
def test_a_malformed_series_says_what_is_wrong(series, message):
    with pytest.raises(SymbolicError, match=message):
        adr_problem(_poisson_with(
            {"region": "all_boundaries", "variable": "u", "type": "neumann",
             "series": series}))




def test_a_robin_condition_becomes_its_three_parts():
    """alpha u + beta (q.n) = g, with g given like any right-hand side."""
    problem = adr_problem(_poisson_with(
        {"region": "all_boundaries", "variable": "u", "type": "robin",
         "alpha": "1.0", "beta": "2.0", "formula": "y**2 - 1"}))
    assert problem["robin"] == [{"component": "u", "region": "all_boundaries",
                                 "where": None, "alpha": "1.0", "beta": "2.0",
                                 "g": "y**2 - 1"}]
    assert problem["dirichlet"] == [] and problem["neumann"] == []


def test_a_robin_condition_needs_both_coefficients():
    with pytest.raises(SymbolicError, match="give alpha and beta"):
        adr_problem(_poisson_with(
            {"region": "all_boundaries", "variable": "u", "type": "robin",
             "alpha": "1.0", "value": 0.0}))


def test_a_robin_condition_needs_its_right_hand_side():
    with pytest.raises(SymbolicError, match="give g"):
        adr_problem(_poisson_with(
            {"region": "all_boundaries", "variable": "u", "type": "robin",
             "alpha": "1.0", "beta": "1.0"}))


def test_the_heat_example_is_transient_with_an_initial_condition():
    problem = adr_problem(load_ymf(EXAMPLES / "heat_equation.ymf"))
    # the time interval comes from the domain, "t ∈ (0, 50]"
    assert problem["time"] == [0.0, 50.0]
    assert problem["has_mass"] is True
    assert problem["initial"]["T"] == "10*numpy.sin(numpy.pi*x)*numpy.sin(numpy.pi*y) + 300"
    # κ came from the coefficients, numerically
    assert problem["adr"]["equations"][0]["diffusion"]["T"]["a"] == ["0.001", "0.001"]


def test_a_composed_navier_stokes_problem_with_periodic_ends():
    problem = adr_problem(load_ymf(EXAMPLES / "navier_stokes" / "plane_poiseuille_re100.ymf"))
    assert problem["unknowns"] == {"v": ["v_0", "v_1"], "p": ["p"]}
    assert {(e["component"], tuple(e["axes"])) for e in problem["periodic"]} == {
        ("v_0", (0,)), ("v_1", (0,)), ("p", (0,))}
    # the forcing (G, 0) with G overridden to 0.08 is a reaction -0.08
    assert problem["adr"]["equations"][0]["reaction"]["r"] == "-0.08"
    # μ overridden to 0.01, so the exact centreline speed stays 1
    assert problem["exact"]["v_0"] == "4.0*y*(1.0 - y)"
    walls = [e for e in problem["dirichlet"] if e["region"] == "walls"]
    assert walls[0]["where"] == "(abs(y) <= tol) or (abs(y - 1.0) <= tol)"   # y = H


def test_coefficient_formulas_may_use_each_other_in_any_order():
    from ymf.symbolic.problem import _coefficients
    values = _coefficients({"f": {"formula": "(G*H, 0)"}, "G": {"value": "2.0"},
                            "H": {"formula": "G + 1"}}, 2)
    assert list(values["f"]) == [6.0, 0]


def test_a_periodic_region_must_hold_both_faces():
    from ymf.symbolic.problem import _periodic_axes
    space = Space(2, {"u": 0})
    with pytest.raises(SymbolicError, match="only one face of axis x"):
        _periodic_axes("x = 0", [0, 0], [4, 1], space, "ends")
    assert _periodic_axes("x = 0 or x = 4", [0, 0], [4, 1], space, "ends") == [0]


def test_dx_of_dx_is_a_divergence():
    d = classify("dx(dx(u)) = 0", {"u": 0})["equations"][0]["diffusion"]["u"]
    assert (d["rowptr"], d["colind"], d["a"]) == ([0, 1, 1], [0], ["-1"])


def test_tensor_conventions():
    space = Space(2, {"v": 1})
    g = parse_expression("grad(v)", space)
    # grad(v)[i, j] = d v_j / d x_i
    assert str(g[0, 1]) == "D(v_1, 0)"
    # dot(v, grad(v)) is (v . nabla) v
    conv = parse_expression("dot(v, grad(v))", space)
    assert str(conv[0]) == "v_0*D(v_0, 0) + v_1*D(v_0, 1)"


def test_tuples_are_vectors_inside_expressions():
    space = Space(2, {"u": 0})
    value = parse_expression("u**2*(1, 0) + (x, 2*y)", space)
    assert list(value) == [sympy.Symbol("u", real=True)**2 + space.x[0], 2 * space.x[1]]


L_SHAPE = {"kind": "piecewise_linear", "dimension": 2,
           "vertices": [[0.0, 0.0], [2.0, 0.0], [2.0, 1.0], [1.0, 1.0],
                        [1.0, 2.0], [0.0, 2.0]],
           "segments": [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]]}


def test_a_piecewise_linear_region_is_carried_with_its_bounding_box():
    doc = load_ymf(EXAMPLES / "poisson.ymf")
    doc["model"]["strong_form"]["domain"] = L_SHAPE
    problem = adr_problem(doc)
    assert problem["domain"]["lower"] == [0.0, 0.0]
    assert problem["domain"]["upper"] == [2.0, 2.0]
    assert problem["domain"]["region"] == L_SHAPE
    assert problem["dirichlet"][0]["where"] is None      # the whole boundary
    yaml.safe_dump(problem)                               # still plain data


def test_periodic_faces_need_a_box():
    doc = _poisson_with({"region": "all_boundaries", "variable": "u", "type": "periodic"})
    doc["model"]["strong_form"]["domain"] = L_SHAPE
    with pytest.raises(SymbolicError, match="periodic faces are a box's"):
        adr_problem(doc)


def test_strong_form_equations_are_what_adr_problem_classifies():
    """The Space and residuals, in sympy, for a tool that renders a spec."""
    from pathlib import Path
    from ymf.compose import load_composed
    from ymf.symbolic import adr_problem, strong_form_equations, to_notation
    doc = load_composed(Path(__file__).parent.parent / "examples" / "kovasznay_flow.ymf")
    space, residuals = strong_form_equations(doc)
    assert space.components == ["v_0", "v_1", "p"]
    assert len(residuals) == 3
    assert [e["component"] for e in adr_problem(doc)["adr"]["equations"]] == space.components
    assert to_notation(residuals[2]) == "∂v_0/∂x + ∂v_1/∂y"
