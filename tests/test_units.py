import pytest

from ydmf.units import (
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
