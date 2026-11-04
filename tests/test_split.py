"""ymf.symbolic.split: each sub-problem of a split discretization, classified.

The two-phase dam break (tests/fixtures/two_phase) is not one ADR system,
but each of its five sub-problems is once what it reads from the others is
data: flow is variable-density Navier-Stokes with ρ and μ read from VOF's
carrier, VOF and the level set are transport in a lagged velocity,
redistancing is an eikonal Hamiltonian, mass correction a regularized
reaction-diffusion.
"""

import copy
import math
from pathlib import Path

import pytest

sympy = pytest.importorskip("sympy")
pytest.importorskip("strictyaml")

from ymf.schema import load_ymf  # noqa: E402
from ymf.symbolic import Space, SymbolicError, adr_form, parse_equation  # noqa: E402
from ymf.symbolic.split import smoothed_heaviside, split_problem  # noqa: E402

TWO_PHASE = Path(__file__).parent / "fixtures" / "two_phase"


@pytest.fixture(scope="module")
def dambreak():
    return load_ymf(TWO_PHASE / "dambreak.ymf")


@pytest.fixture(scope="module")
def split(dambreak):
    return split_problem(dambreak, "level_set_five")


def _sub(split, name):
    return next(s for s in split["subproblems"] if s["name"] == name)


def _equation(sub, component):
    return next(e for e in sub["adr"]["equations"] if e["component"] == component)


def _evaluate(code, **names):
    numpy = pytest.importorskip("numpy")       # not in ymf[symbolic] alone
    return eval(code, {"__builtins__": {}, "numpy": numpy}, names)  # noqa: S307


# --- the setting ---------------------------------------------------------------


def test_the_setting(split):
    assert split["domain"] == {"lower": [0.0, 0.0], "upper": [3.22, 1.8]}
    assert split["time"] == [0.0, 3.0]
    assert split["bound"] == {"h": "element_size", "ε": "1.5*h"}
    assert set(split["initial"]) == {"v_0", "v_1", "p", "φ"}


def test_each_sub_problem_has_its_own_unknowns(split):
    assert [(s["name"], s["unknowns"]) for s in split["subproblems"]] == [
        ("flow", {"v": ["v_0", "v_1"], "p": ["p"]}),
        ("vof", {"χh": ["χh"]}),
        ("ncls", {"φ": ["φ"]}),
        ("rdls", {"φ": ["φ"]}),
        ("mcorr", {"φc": ["φc"]}),
    ]


# --- the five sub-problems -------------------------------------------------------


def test_flow_reads_density_and_viscosity_through_the_vof_carrier(split):
    flow = _sub(split, "flow")
    assert flow["reads"]["χ"] == {"at": "previous_step", "via": "vof", "carried_by": "χh",
                                  "components": ["χh"]}
    assert flow["coefficients"]["ρ"]["code"] == "998.2 - 996.995*χh"
    assert flow["coefficients"]["ρ"]["reads"] == ["χh"]
    momentum = _equation(flow, "v_0")
    assert momentum["mass"]["m"] == "v_0*(998.2 - 996.995*χh)"
    assert set(momentum["diffusion"]) == {"v_0", "v_1"}          # 2μ sym(∇v)
    assert flow["solver"]["max_iterations"] == 100                 # refined
    assert flow["proteus"] == {"model": "RANS2P", "useVF": "1.0"}  # a hint, carried as text
    assert "regularized" not in flow           # χ is read through χh, not evaluated


def test_the_level_set_is_a_hamiltonian_in_the_lagged_velocity(split):
    ncls = _sub(split, "ncls")
    assert ncls["reads"] == {"v": {"at": "current_step", "via": "flow",
                                   "components": ["v_0", "v_1"]}}
    h = _equation(ncls, "φ")["hamiltonian"]
    assert h["H"] == "grad_φ_0*v_0 + grad_φ_1*v_1"
    assert h["dH"] == {"φ": ["v_0", "v_1"]}


def test_vof_is_conservative_transport(split):
    vof = _equation(_sub(split, "vof"), "χh")
    assert vof["advection"]["f"] == ["χh*v_0", "χh*v_1"]


def test_vof_starts_from_the_smoothed_level_set(split):
    vof = _sub(split, "vof")
    assert vof["initial"]["reads"] == ["φ"]
    assert vof["regularized"] == {"heaviside": {"by": "smoothed_heaviside", "width": "1.5*h"}}
    assert _evaluate(vof["initial"]["χh"], φ=0.0, h=0.1) == 0.5
    assert _evaluate(vof["initial"]["χh"], φ=1.0, h=0.1) == 1.0


def test_redistancing_reads_a_copy_of_the_level_set_and_smooths_its_sign(split):
    rdls = _sub(split, "rdls")
    assert rdls["reads"]["φs"]["of"] == "φ"
    assert rdls["regularized"] == {"sign": {"by": "smoothed_sign", "width": "0.75*h"}}
    eq = _equation(rdls, "φ")
    assert eq["hamiltonian"]["depends_on"] == {"φ": "nonlinear"}  # |∇φ|
    # the reaction is -S(φs): 1 - 2 H_ε(φs)
    assert _evaluate(eq["reaction"]["r"], φs=1.0, h=0.1) == -1.0
    assert _evaluate(eq["reaction"]["r"], φs=-1.0, h=0.1) == 1.0


def test_mass_correction_is_a_regularized_reaction_diffusion(split):
    mcorr = _sub(split, "mcorr")
    eq = _equation(mcorr, "φc")
    assert eq["reaction"]["depends_on"] == {"φc": "nonlinear"}
    assert eq["diffusion"]["φc"]["a"] == ["10*h", "10*h"]
    # at φ + φc = 0 the reaction is H_ε(0) - χh = 1/2 - χh, its derivative 1/ε
    assert _evaluate(eq["reaction"]["r"], φ=0.0, φc=0.0, χh=0.25, h=0.1) == 0.25
    dr = eq["reaction"]["dr"]["φc"]
    assert math.isclose(_evaluate(dr, φ=0.0, φc=0.0, χh=0.25, h=0.1), 1 / 0.15)


# --- the regularization is Proteus's ----------------------------------------------


def _proteus_H(eps, phi):           # equivalent_polynomials.h, gf.H
    if phi > eps:
        return 1.0
    if phi < -eps:
        return 0.0
    if phi == 0.0:
        return 0.5
    return 0.5 * (1.0 + phi / eps + math.sin(math.pi * phi / eps) / math.pi)


def _proteus_D(eps, phi):           # gf.D
    if phi > eps or phi < -eps:
        return 0.0
    return 0.5 * (1.0 + math.cos(math.pi * phi / eps)) / eps


def test_the_smoothed_heaviside_and_its_derivative_are_proteus_gf_H_and_gf_D():
    from ymf.symbolic.adr import to_code
    phi, eps = sympy.symbols("phi eps", real=True)
    H = to_code(smoothed_heaviside(phi, eps))
    D = to_code(sympy.diff(smoothed_heaviside(phi, eps), phi))
    for e in (0.01, 0.0483, 0.3):
        for k in range(-40, 41):
            p = 2 * e * k / 40
            assert math.isclose(_evaluate(H, phi=p, eps=e), _proteus_H(e, p), abs_tol=1e-15)
            assert math.isclose(_evaluate(D, phi=p, eps=e), _proteus_D(e, p),
                                rel_tol=1e-13, abs_tol=1e-13)


# --- what is refused -----------------------------------------------------------------


def _without(dambreak, sub, key):
    doc = copy.deepcopy(dambreak)
    entry = doc["solution_paths"]["discretizations"][0]
    next(s for s in entry["subproblems"] if s["name"] == sub).pop(key)
    return doc


def test_a_field_read_without_saying_when_is_refused(dambreak):
    with pytest.raises(SymbolicError, match="sub-problem ncls.*declared in its lags"):
        split_problem(_without(dambreak, "ncls", "lags"), "level_set_five")


def test_an_unregularized_sharp_function_of_an_unknown_is_refused(dambreak):
    doc = copy.deepcopy(dambreak)
    doc["solution_paths"]["discretizations"][0].pop("regularize")
    with pytest.raises(SymbolicError, match="sub-problem mcorr.*regularizes it"):
        split_problem(doc, "level_set_five")


def test_a_regularization_ymf_does_not_know_is_refused(dambreak):
    doc = copy.deepcopy(dambreak)
    doc["solution_paths"]["discretizations"][0]["regularize"]["heaviside"]["by"] = "tanh"
    with pytest.raises(SymbolicError, match="regularized by 'tanh'"):
        split_problem(doc, "level_set_five")


def test_only_split_discretizations_and_real_names(dambreak):
    with pytest.raises(SymbolicError, match="no discretization named 'nope'"):
        split_problem(dambreak, "nope")
    doc = copy.deepcopy(dambreak)
    doc["solution_paths"]["discretizations"][0].pop("subproblems")
    with pytest.raises(SymbolicError, match="not split into sub-problems"):
        split_problem(doc, "level_set_five")


def test_a_sub_problem_that_only_spins_up_is_carried_unclassified(dambreak):
    doc = copy.deepcopy(dambreak)
    doc["solution_paths"]["discretizations"][0]["subproblems"].append(
        {"name": "warmup", "spin_up": True})
    warmup = _sub(split_problem(doc, "level_set_five"), "warmup")
    assert warmup["classified"] is False and warmup["why"] == "it spins the run up"


# --- known fields in the language ----------------------------------------------------


def test_a_known_field_is_differentiated_as_a_field():
    space = Space(2, {"c": 0}, {}, fields={"v": 1})
    reaction = adr_form(parse_equation("∂c/∂t = ∇·v", space), space)["equations"][0]["reaction"]
    assert reaction["r"] == "-grad_v_0_0 - grad_v_1_1"


def test_the_rate_of_a_known_field_has_no_name():
    space = Space(2, {"c": 0}, {}, fields={"v": 1})
    with pytest.raises(SymbolicError, match="time derivative of the known field v_0"):
        adr_form(parse_equation("∂c/∂t = ∂v_0/∂t", space), space)


def test_a_name_is_an_unknown_or_a_known_field_not_both():
    with pytest.raises(SymbolicError, match="both an unknown and a known field"):
        Space(2, {"v": 1}, {}, fields={"v": 1})


# --- boundary conditions, per sub-problem -----------------------------------------------


def test_each_sub_problem_imposes_the_conditions_on_what_it_solves(split):
    flow = [(c["type"], c["variable"], c["region"], c.get("set"))
            for c in _sub(split, "flow")["boundary_conditions"]]
    assert flow == [("slip", "v", "walls", "free_slip"),
                    ("dirichlet", "p", "top", "atmosphere"),
                    ("open", "v", "top", "atmosphere")]
    slip = _sub(split, "flow")["boundary_conditions"][0]
    assert slip["components"] == ["v_0", "v_1"] and slip["value"] == "0.0"
    assert _sub(split, "ncls")["boundary_conditions"] == []


def test_a_condition_on_a_coefficient_goes_to_its_carrier(split):
    (inflow,) = _sub(split, "vof")["boundary_conditions"]
    assert inflow == {"type": "inflow", "variable": "χ", "region": "top",
                      "where": "(abs(y - 1.8) <= tol)", "carried_by": "χh",
                      "set": "atmosphere", "components": ["χh"], "value": "1.0"}


def test_a_condition_reaches_every_sub_problem_that_solves_for_its_variable(dambreak):
    doc = copy.deepcopy(dambreak)
    doc["model"]["strong_form"]["boundary_conditions"] = [
        {"region": "top", "variable": "φ", "type": "dirichlet", "value": 1.0},
        {"region": "top", "variable": "v", "type": "inflow", "formula": "(0, -t)"},
        {"region": "walls", "variable": "p", "type": "characteristic",
         "series": [[0.0, 0.0], [1.0, 2.0]]}]
    split = split_problem(doc, "level_set_five")
    on_phi = [s["name"] for s in split["subproblems"]
              if any(c["variable"] == "φ" for c in s["boundary_conditions"])]
    assert on_phi == ["ncls", "rdls"]                 # not mcorr, which only reads φ
    flow = {(c["variable"], c["region"]): c for c in _sub(split, "flow")["boundary_conditions"]}
    assert flow[("v", "top")]["value"] == {"v_0": "0", "v_1": "-t"}   # replaces open
    assert flow[("p", "walls")]["series"] == [[0.0, 0.0], [1.0, 2.0]]
