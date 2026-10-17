import math
import pytest

# The YMF core deliberately depends on pyyaml alone, so a solver writing
# archives need not carry the validation stack. This whole module is about
# ymf.units, which is the pint-backed half of that stack.
pytest.importorskip("pint")

from ymf.units import (
    UnitCheckResult,
    UnitParseError,
    characteristic_scale_values,
    check,
    compute_dimensionless_numbers,
    non_dimensionalize,
    normalize_unit_string,
    to_pint_unit,
    units_compatible,
)


# ---------------------------------------------------------------------------
# normalize_unit_string / to_pint_unit
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "raw,expected_normalized",
    [
        ("m/s", "m/s"),
        ("m2/s", "m**2/s"),
        ("kg/m3", "kg/m**3"),
        ("kg/(m*s)", "kg/(m*s)"),
        ("kg/(m·s)", "kg/(m*s)"),
        ("m^2", "m**2"),
        ("m²", "m**2"),
    ],
)
def test_normalize_unit_string(raw, expected_normalized):
    assert normalize_unit_string(raw) == expected_normalized


def test_to_pint_unit_none_and_dimensionless():
    assert to_pint_unit(None).dimensionless
    assert to_pint_unit("dimensionless").dimensionless
    assert to_pint_unit("").dimensionless


def test_to_pint_unit_bad_string_raises():
    with pytest.raises(UnitParseError):
        to_pint_unit("not_a_unit_at_all_xyz")


def test_units_compatible():
    assert units_compatible("m/s", "km/hour")
    assert not units_compatible("m/s", "Pa")


# ---------------------------------------------------------------------------
# check()
# ---------------------------------------------------------------------------


def _base_doc(**overrides):
    doc = {
        "model": {
            "name": "test",
            "strong_form": {
                "unknowns": [
                    {"name": "u", "units": "m/s"},
                    {"name": "p", "units": "Pa"},
                ],
                "boundary_conditions": [
                    {"region": "inlet", "variable": "u", "type": "dirichlet",
                     "value": 5.0, "units": "m/s"},
                ],
                "initial_conditions": [
                    {"field": "u", "type": "constant", "value": 0.0, "units": "m/s"},
                ],
            },
        },
        "solution_paths": {"analytical": [], "discretizations": []},
    }
    doc["model"]["strong_form"].update(overrides)
    return doc


def test_check_no_units_at_all_is_ok():
    doc = {
        "model": {"name": "bare", "strong_form": {"unknowns": ["u", "p"]}},
        "solution_paths": {"analytical": [], "discretizations": []},
    }
    result = check(doc)
    assert isinstance(result, UnitCheckResult)
    assert result.ok
    assert bool(result) is True
    assert len(result) == 0


def test_check_consistent_units_ok():
    result = check(_base_doc())
    assert result.ok, [str(i) for i in result]


def test_check_bc_unit_mismatch_flagged():
    doc = _base_doc(
        boundary_conditions=[
            {"region": "inlet", "variable": "u", "type": "dirichlet",
             "value": 5.0, "units": "Pa"},  # wrong: u is m/s
        ]
    )
    result = check(doc)
    assert not result.ok
    assert any("incompatible" in str(i) for i in result)


def test_check_ic_unit_mismatch_flagged():
    doc = _base_doc(
        initial_conditions=[
            {"field": "u", "type": "constant", "value": 0.0, "units": "kg"},
        ]
    )
    result = check(doc)
    assert not result.ok


def test_check_bad_unit_string_on_unknown_flagged():
    doc = _base_doc(unknowns=[{"name": "u", "units": "totally_bogus_unit"}])
    result = check(doc)
    assert not result.ok


# ---------------------------------------------------------------------------
# characteristic_scale_values / compute_dimensionless_numbers
# ---------------------------------------------------------------------------


def test_characteristic_scale_values_direct_and_derived():
    problem = {
        "characteristic_scales": {
            "length": {"value": 1.0, "units": "m"},
            "velocity": {"value": 5.0, "units": "m/s"},
            "time": {"derived": True, "formula": "L_char / U_char"},
        }
    }
    scales = characteristic_scale_values(problem)
    assert scales["length"] == 1.0
    assert scales["velocity"] == 5.0
    assert scales["time"] == pytest.approx(0.2)


def test_compute_dimensionless_numbers_reynolds():
    problem = {
        "characteristic_scales": {
            "length": {"value": 1.0, "units": "m"},
            "velocity": {"value": 5.0, "units": "m/s"},
            "density": {"value": 1.225, "units": "kg/m3"},
            "viscosity": {"value": 1.81e-5, "units": "kg/(m*s)"},
        },
        "dimensionless_numbers": {
            "Reynolds": {"formula": "rho * U * L / mu"},
        },
    }
    numbers = compute_dimensionless_numbers(problem)
    assert numbers["Reynolds"] == pytest.approx((1.225 * 5.0 * 1.0) / 1.81e-5)


def test_compute_dimensionless_numbers_missing_scale_omitted():
    problem = {
        "characteristic_scales": {"length": {"value": 1.0, "units": "m"}},
        "dimensionless_numbers": {
            "Reynolds": {"formula": "rho * U * L / mu"},  # rho/U/mu missing
        },
    }
    numbers = compute_dimensionless_numbers(problem)
    assert "Reynolds" not in numbers


def _heat_doc_with_numbers(coefficients_yaml):
    """examples/heat_equation.ymf with a velocity scale, the given
    coefficients block, and coefficient-referencing dimensionless numbers."""
    from pathlib import Path

    text = (Path(__file__).parent.parent / "examples" / "heat_equation.ymf").read_text()
    scales = '    time: {derived: true, formula: "1/(2π²κ)", units: s}\n'
    coefficients = "    coefficients:\n      κ: {value: 1.0e-3, units: m2/s}\n"
    assert scales in text and coefficients in text
    text = text.replace(
        scales,
        scales
        + "    velocity: {value: 0.5, units: m/s}\n"
        + "\n  dimensionless_numbers:\n"
        + '    Peclet: {formula: "U * L / κ"}\n'
        + '    Fourier: {formula: "κ * time / L**2"}\n',
    )
    return text.replace(coefficients, "    coefficients:\n" + coefficients_yaml)


@pytest.mark.parametrize(
    "coefficients_yaml",
    [
        "      κ: {value: 1.0e-3, units: m2/s}\n",
        "      κ: 1.0e-3\n",
        '      κ: "1.0e-3"\n',
    ],
    ids=["value-map", "bare-float", "quoted-float"],
)
def test_compute_dimensionless_numbers_uses_validated_coefficients(coefficients_yaml):
    # Regression: coefficients are MapPattern(Str(), Any()) in the schema, so
    # strictyaml returns '1.0e-3' (a string) and formulas referencing kappa
    # used to be silently dropped.
    from ymf import validate_ymf

    doc = validate_ymf(_heat_doc_with_numbers(coefficients_yaml))
    numbers = compute_dimensionless_numbers(doc["model"])
    assert numbers["Peclet"] == pytest.approx(0.5 * 1.0 / 1.0e-3)
    # the time scale is derived, 1/(2π²κ), so κ·time/L² is 1/(2π²) exactly
    assert numbers["Fourier"] == pytest.approx(1.0 / (2.0 * math.pi ** 2))


def test_an_expression_coefficient_resolves_but_a_non_number_still_does_not():
    """"π²" is a number once evaluated; "nan" and True are not numbers at all.

    This used to drop all three, so a dimensionless number built on an
    expression-valued coefficient went silently missing while adr_problem
    resolved the same coefficient without difficulty.
    """
    problem = {
        "characteristic_scales": {"length": {"value": 1.0, "units": "m"}},
        "strong_form": {"coefficients": {"Q": "π²", "k": "nan", "flag": True}},
        "dimensionless_numbers": {
            "A": {"formula": "Q * L"},
            "B": {"formula": "k * L"},
            "C": {"formula": "flag * L"},
        },
    }
    assert compute_dimensionless_numbers(problem) == {
        "A": pytest.approx(math.pi ** 2)}


def test_units_and_the_solver_agree_on_a_coefficient():
    """One spec, one answer -- the reason sympy is in the frontend extra."""
    from ymf import load_ymf
    from ymf.symbolic.problem import coefficient_numbers
    from ymf.units import _coefficient_values

    pytest.importorskip("sympy")
    problem = load_ymf("examples/kovasznay_flow.ymf")["model"]
    by_units = _coefficient_values(problem)
    by_solver = coefficient_numbers(problem["strong_form"]["coefficients"])
    assert by_units == by_solver
    # λ is written as an expression over ρ and μ; both sides evaluate it
    assert by_units["λ"] == pytest.approx(-0.9637405441957657)


# ---------------------------------------------------------------------------
# non_dimensionalize()
# ---------------------------------------------------------------------------


def test_non_dimensionalize_produces_substitutions_for_velocity_and_pressure():
    doc = {
        "model": {
            "name": "ns",
            "characteristic_scales": {
                "length": {"value": 1.0, "units": "m"},
                "velocity": {"value": 10.0, "units": "m/s"},
                "density": {"value": 1.225, "units": "kg/m3"},
                "viscosity": {"value": 1.5e-5, "units": "kg/(m*s)"},
            },
            "dimensionless_numbers": {
                "Reynolds": {"formula": "rho * U * L / mu"},
            },
            "strong_form": {
                "unknowns": [
                    {"name": "u", "units": "m/s"},
                    {"name": "p", "units": "Pa"},
                ],
            },
        }
    }
    result = non_dimensionalize(doc)
    assert result["substitutions"]["u*"] == "u / U_char"
    assert result["substitutions"]["p*"] == "p / (rho_char * U_char**2)"
    assert result["dimensionless_numbers"]["Reynolds"] == pytest.approx(
        (1.225 * 10.0 * 1.0) / 1.5e-5
    )


def test_non_dimensionalize_no_units_produces_empty_substitutions():
    doc = {"model": {"name": "bare", "strong_form": {"unknowns": ["u", "p"]}}}
    result = non_dimensionalize(doc)
    assert result["substitutions"] == {}


# --- scales derived from coefficients ---------------------------------------


def _flow(scales, coefficients):
    return {"model": {
        "characteristic_scales": scales,
        "dimensionless_numbers": {"Re": {"formula": "rho * U * L / mu"}},
        "strong_form": {"unknowns": [], "coefficients": coefficients},
    }}


def test_a_scale_can_be_derived_from_coefficients():
    from ymf.units import characteristic_scale_values, compute_dimensionless_numbers

    # strings, as strictyaml returns coefficient values
    doc = _flow({"length": {"value": 1.0},
                 "velocity": {"derived": True, "formula": "G * L_char**2 / (8 * mu)"},
                 "density": {"derived": True, "formula": "rho"},
                 "viscosity": {"derived": True, "formula": "mu"}},
                {"rho": {"value": "1.0"}, "mu": {"value": "0.01"}, "G": {"value": "0.08"}})
    scales = characteristic_scale_values(doc["model"])
    assert scales["viscosity"] == 0.01
    assert scales["velocity"] == pytest.approx(1.0)
    assert compute_dimensionless_numbers(doc["model"])["Re"] == pytest.approx(100.0)


def test_a_derived_scale_raises_no_duplication_issue():
    doc = _flow({"viscosity": {"derived": True, "formula": "mu"}}, {"mu": {"value": "0.01"}})
    assert list(check(doc)) == []


def test_a_scale_restating_a_coefficient_with_a_different_value_is_an_error():
    doc = _flow({"viscosity": {"value": 1.0}}, {"mu": {"value": "0.01"}})
    result = check(doc)
    assert not result
    (issue,) = list(result)
    assert issue.severity == "error"
    assert issue.location == "characteristic_scales[viscosity]"
    assert 'formula: "mu"' in issue.message


def test_a_scale_restating_a_coefficient_with_the_same_value_is_a_warning():
    doc = _flow({"density": {"value": 1.0}}, {"rho": 1.0})
    result = check(doc)
    assert result
    assert [i.severity for i in result] == ["warning"]


# ---------------------------------------------------------------------------
# one value written twice, found by value rather than by name
# ---------------------------------------------------------------------------
#
# Phase 4 of the semantics campaign. A scale stated as a number beside a
# coefficient for the same quantity is one value written twice, and the
# copies drift when someone changes one. The check matched the conventional
# names only -- rho/ρ for density, mu/μ for viscosity -- so a coefficient
# called nu or eta, or any scale outside the fixed six, was invisible to it.
#
# Matching on the value needs something to tell a duplicate from a
# coincidence, because in a non-dimensionalized problem half the quantities
# are 1.0. Dimension is that something.


def _problem(coefficients, scales):
    return {"model": {
        "name": "p",
        "strong_form": {
            "provenance": "llm_derived", "equation_formulation": "f",
            "unknowns": ["u"], "domain": "Ω",
            "boundary_regions": [{"name": "a", "where": "∂Ω"}],
            "coefficients": coefficients},
        "characteristic_scales": scales}}


def _scale_issues(coefficients, scales):
    from ymf.units import check

    return [i for i in check(_problem(coefficients, scales)).issues
            if "characteristic_scales" in i.location]


def test_a_coefficient_under_an_unconventional_name_is_now_found():
    """nu, not mu -- the case the name matching could never see."""
    (issue,) = _scale_issues(
        {"nu": {"value": 0.025, "units": "Pa*s"}},
        {"viscosity": {"value": 0.025, "units": "Pa*s"}})
    assert issue.severity == "warning"
    assert "equals coefficients[nu]" in issue.message
    assert 'formula: "nu"' in issue.message


def test_a_scale_outside_the_conventional_six_is_covered():
    """characteristic_scales takes any name now, so the check has to as well."""
    (issue,) = _scale_issues(
        {"kappa": {"value": 1.0e-3, "units": "m2/s"}},
        {"diffusivity": {"value": 1.0e-3, "units": "m2/s"}})
    assert "equals coefficients[kappa]" in issue.message


def test_equal_numbers_of_different_kinds_are_a_coincidence():
    """1.0 kg/m3 and 1.0 m are not one value written twice.

    Without this, every non-dimensionalized problem would report a handful
    of duplicates that are nothing of the sort.
    """
    assert _scale_issues(
        {"rho_": {"value": 1.0, "units": "kg/m3"}},
        {"length": {"value": 1.0, "units": "m"}}) == []


def test_nothing_is_said_when_either_side_states_no_units():
    """With no dimension there is no way to tell duplication from accident."""
    assert _scale_issues(
        {"nu": {"value": 0.025}},
        {"viscosity": {"value": 0.025, "units": "Pa*s"}}) == []


def test_the_conventional_pairing_still_reports_a_disagreement():
    """The name match is what catches copies that have already drifted.

    Value matching cannot: once the copies differ, the values no longer
    match, which is exactly when it matters most.
    """
    (issue,) = _scale_issues(
        {"mu": {"value": 0.1, "units": "Pa*s"}},
        {"viscosity": {"value": 0.025, "units": "Pa*s"}})
    assert issue.severity == "error"
    assert "differs from" in issue.message


def test_a_pair_is_reported_once_not_twice():
    """mu matches by name and by value; it is one finding."""
    issues = _scale_issues(
        {"mu": {"value": 0.025, "units": "Pa*s"}},
        {"viscosity": {"value": 0.025, "units": "Pa*s"}})
    assert len(issues) == 1


# ---------------------------------------------------------------------------
# boundary conditions: the units of g against the type of condition
# ---------------------------------------------------------------------------


def _conduction(bc, unknown_units="K", domain_units="m", k_units="W/(m*K)",
                h_units="W/(m**2*K)"):
    """Steady conduction, -div(k grad T) = 0, with one boundary condition.

    The diffusive flux is q = -k grad T, so [q] = [k][T]/[L], W/m2 here.
    """
    domain = {"kind": "box", "lower": [0.0, 0.0], "upper": [1.0, 1.0]}
    if domain_units:
        domain["units"] = domain_units
    k = {"value": 2.0}
    if k_units:
        k["units"] = k_units
    h = {"value": 10.0}
    if h_units:
        h["units"] = h_units
    unknown = {"name": "T"}
    if unknown_units:
        unknown["units"] = unknown_units
    return {
        "model": {"name": "conduction", "strong_form": {
            "unknowns": [unknown],
            "equations": ["-∇·(k∇T) = 0  in Ω"],
            "domain": domain,
            "boundary_regions": [{"name": "wall", "where": "x = 0"}],
            "boundary_conditions": [dict({"region": "wall", "variable": "T"}, **bc)],
            "coefficients": {"k": k, "h": h},
        }},
        "solution_paths": {"analytical": [], "discretizations": []},
    }


def _robin(**kw):
    """h T + k q.n = g: [h][T] = W/m2 and [k][q] = W2/(m3 K) -- not equal, so
    beta is k/k, or 1, in the consistent form below."""
    return _conduction(dict({"type": "robin", "alpha": "h", "beta": "1",
                             "value": 3000.0, "units": "W/m**2"}, **kw.pop("bc", {})), **kw)


def _errors(result):
    return [str(i) for i in result if i.severity == "error"]


def test_a_neumann_g_is_a_flux_not_a_value_of_the_unknown():
    """The check used to compare every condition's units with the unknown's,
    which is right for dirichlet only: a neumann g in K would have passed
    and one in W/m2 would have failed."""
    pytest.importorskip("ymf.symbolic")
    ok = check(_conduction({"type": "neumann", "value": 5.0, "units": "W/m**2"}))
    assert ok.ok, _errors(ok)
    wrong = check(_conduction({"type": "neumann", "value": 5.0, "units": "K"}))
    assert not wrong.ok
    assert any("flux q.n of T" in e for e in _errors(wrong))


def test_a_dirichlet_g_is_still_a_value_of_the_unknown():
    assert check(_conduction({"type": "dirichlet", "value": 300.0, "units": "degC"})).ok
    wrong = check(_conduction({"type": "dirichlet", "value": 300.0, "units": "W/m**2"}))
    assert any("incompatible" in e for e in _errors(wrong))


def test_a_consistent_robin_condition_passes():
    pytest.importorskip("ymf.symbolic")
    result = check(_robin())
    assert result.ok and len(result) == 0, [str(i) for i in result]


def test_robin_alpha_against_beta_needs_no_units_for_t_or_g():
    """[alpha][u] = [beta][q.n] holds or fails whatever u's units are: u is a
    factor of both sides. Here beta = 1 makes beta q.n W/m2 while alpha T is
    h's W/(m2 K) times T, so they agree only if T is K -- undecided without
    T's units. With beta = k/h they agree whatever T is."""
    pytest.importorskip("ymf.symbolic")
    bare = dict(unknown_units=None, bc={"units": None})
    undecided = check(_robin(**bare), missing_units="error")
    assert not undecided.ok
    assert all("not checked" in e for e in _errors(undecided))
    # alpha = h, beta = k (W/(m K)): k q.n is W2/(m3 K), h T is W T/(m2 K)
    mismatched = check(_robin(unknown_units=None,
                              bc={"beta": "k", "units": None}), missing_units="error")
    assert any("alpha u is" in e and "beta q.n is" in e for e in _errors(mismatched))
    # alpha = h/k (1/m), beta = 1: (1/m) T against q/... -> [T]/m = [k][T]/[L]/[k]
    consistent = check(_robin(unknown_units=None,
                              bc={"alpha": "h/k", "beta": "1/k", "units": None}),
                       missing_units="error")
    errors = _errors(consistent)
    # alpha u ~ beta q.n is decided without T; g's relation to them is not
    assert errors and all("units of T" in e or "units of g" in e for e in errors)
    assert not any("alpha u is" in e for e in errors)


def test_a_robin_g_that_disagrees_is_an_error():
    pytest.importorskip("ymf.symbolic")
    result = check(_robin(bc={"units": "K"}))
    assert any(e.startswith("[error]") and "g is" in e for e in _errors(result))


def test_missing_units_are_skipped_by_default():
    """Most specs declare no units; checking them must say nothing unless asked."""
    pytest.importorskip("ymf.symbolic")
    bare = _conduction({"type": "neumann", "value": 1.0}, unknown_units=None,
                       domain_units=None, k_units=None, h_units=None)
    assert len(check(bare)) == 0
    warned = check(bare, missing_units="warn")
    assert warned.ok and len(warned) == 1
    assert [i.severity for i in warned] == ["warning"]
    assert "not checked: nothing declares" in str(warned.issues[0])
    errored = check(bare, missing_units="error")
    assert not errored.ok


def test_missing_units_names_what_is_undeclared():
    pytest.importorskip("ymf.symbolic")
    result = check(_conduction({"type": "neumann", "value": 1.0, "units": "W/m**2"},
                               k_units=None), missing_units="warn")
    assert [i.severity for i in result] == ["warning"]
    assert "the units of k" in result.issues[0].message
    assert "the units of T" not in result.issues[0].message


def test_dirichlet_missing_units_follow_the_option_too():
    doc = _conduction({"type": "dirichlet", "value": 300.0})   # g has no units
    assert len(check(doc)) == 0
    assert [i.severity for i in check(doc, missing_units="warn")] == ["warning"]


def test_missing_units_must_be_a_known_policy():
    with pytest.raises(ValueError, match="missing_units"):
        check(_conduction({"type": "dirichlet", "value": 300.0}), missing_units="loud")


def test_checking_twice_reuses_the_stand_in_units():
    """Stand-ins are pint units; defining them again on each check would be
    a redefinition in pint's registry."""
    pytest.importorskip("ymf.symbolic")
    doc = _conduction({"type": "neumann", "value": 1.0}, unknown_units=None)
    first = [str(i) for i in check(doc, missing_units="warn")]
    assert first == [str(i) for i in check(doc, missing_units="warn")]
