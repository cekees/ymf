import math
import pytest

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
        "Problem": {
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
    doc["Problem"]["strong_form"].update(overrides)
    return doc


def test_check_no_units_at_all_is_ok():
    doc = {
        "Problem": {"name": "bare", "strong_form": {"unknowns": ["u", "p"]}},
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
    numbers = compute_dimensionless_numbers(doc["Problem"])
    assert numbers["Peclet"] == pytest.approx(0.5 * 1.0 / 1.0e-3)
    # the time scale is derived, 1/(2π²κ), so κ·time/L² is 1/(2π²) exactly
    assert numbers["Fourier"] == pytest.approx(1.0 / (2.0 * math.pi ** 2))


def test_compute_dimensionless_numbers_ignores_non_numeric_coefficients():
    problem = {
        "characteristic_scales": {"length": {"value": 1.0, "units": "m"}},
        "strong_form": {"coefficients": {"Q": "π²", "k": "nan", "flag": True}},
        "dimensionless_numbers": {
            "A": {"formula": "Q * L"},
            "B": {"formula": "k * L"},
            "C": {"formula": "flag * L"},
        },
    }
    assert compute_dimensionless_numbers(problem) == {}


# ---------------------------------------------------------------------------
# non_dimensionalize()
# ---------------------------------------------------------------------------


def test_non_dimensionalize_produces_substitutions_for_velocity_and_pressure():
    doc = {
        "Problem": {
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
    doc = {"Problem": {"name": "bare", "strong_form": {"unknowns": ["u", "p"]}}}
    result = non_dimensionalize(doc)
    assert result["substitutions"] == {}


# --- scales derived from coefficients ---------------------------------------


def _flow(scales, coefficients):
    return {"Problem": {
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
    scales = characteristic_scale_values(doc["Problem"])
    assert scales["viscosity"] == 0.01
    assert scales["velocity"] == pytest.approx(1.0)
    assert compute_dimensionless_numbers(doc["Problem"])["Re"] == pytest.approx(100.0)


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
