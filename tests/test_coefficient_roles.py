"""Coefficient roles (C42 round 3, decision 20): what a coefficient means,
whatever the model calls it, so a solver finds the density, the gravity or
the bed by its role and not by name.
"""

import copy
from pathlib import Path

import pytest

pytest.importorskip("strictyaml")

from ymf import vocabulary as v  # noqa: E402
from ymf.references import YmfReferenceError, resolve  # noqa: E402
from ymf.schema import load_ymf  # noqa: E402

FIXTURES = Path(__file__).parent / "fixtures"


def test_the_two_phase_family_states_its_roles():
    doc = load_ymf(FIXTURES / "two_phase" / "dambreak.ymf")
    assert v.coefficient_roles(doc) == {"phase_indicator": "χ", "density": "ρ",
                                        "dynamic_viscosity": "μ", "gravity": "g"}
    assert v.coefficient_by_role(doc, "density") == "ρ"
    assert v.coefficient_by_role(doc, "surface_tension") is None


def test_the_shallow_water_family_states_its_roles():
    doc = load_ymf(FIXTURES / "shallow_water" / "lake_at_rest.ymf")
    roles = v.coefficient_roles(doc)
    assert roles == {"bed_friction": "τb", "canopy_drag": "Fd", "gravity": "g",
                     "manning_coefficient": "n", "porosity": "θ", "bed_elevation": "Z",
                     "linear_drag": "α", "quadratic_drag": "β"}
    # the dispatch table's names for them are aliases
    assert v.coefficient_by_role(doc, "mannings") == "n"
    assert v.coefficient_by_role(doc, "bathymetry") == "Z"


def test_a_role_is_the_closed_specs_too():
    from ymf.explicit import explicit
    closed = explicit(load_ymf(FIXTURES / "two_phase" / "dambreak.ymf"))
    assert closed["model"]["strong_form"]["coefficients"]["ρ"]["role"] == "density"


def test_a_role_is_used_once_and_must_be_known():
    doc = load_ymf(FIXTURES / "two_phase" / "dambreak.ymf")
    bad = copy.deepcopy(doc)
    bad["model"]["strong_form"]["coefficients"]["ρw"]["role"] = "density"
    with pytest.raises(YmfReferenceError, match="ρw.role: density is already the role of ρ"):
        resolve(bad)
    bad = copy.deepcopy(doc)
    bad["model"]["strong_form"]["coefficients"]["ρw"]["role"] = "water_density"
    with pytest.raises(YmfReferenceError, match="coefficient role water_density: ymf knows density"):
        resolve(bad)


def test_units_must_be_the_roles():
    pytest.importorskip("pint")
    from ymf.units import check
    doc = load_ymf(FIXTURES / "two_phase" / "dambreak.ymf")
    assert check(doc).ok
    bad = copy.deepcopy(doc)
    bad["model"]["strong_form"]["coefficients"]["μ"]["units"] = "m2/s"
    result = check(bad)
    assert not result.ok
    assert any("units m2/s are not those of its role dynamic_viscosity (Pa*s)" in i.message
               for i in result)
