import pytest

from ymf.normalize import normalize_unknowns


def test_normalize_v01_bare_strings():
    result = normalize_unknowns(["u", "p"])
    assert result == [
        {"name": "u", "units": None, "std_name": None, "provenance": None},
        {"name": "p", "units": None, "std_name": None, "provenance": None},
    ]


def test_normalize_v02_enriched_maps():
    result = normalize_unknowns(
        [
            {"name": "u", "units": "m/s", "std_name": "surface_water_velocity_x"},
            {"name": "p", "units": "Pa"},
        ]
    )
    assert result == [
        {"name": "u", "units": "m/s", "std_name": "surface_water_velocity_x", "provenance": None},
        {"name": "p", "units": "Pa", "std_name": None, "provenance": None},
    ]


def test_normalize_mixed_forms():
    result = normalize_unknowns(["u", {"name": "p", "units": "Pa"}])
    assert result == [
        {"name": "u", "units": None, "std_name": None, "provenance": None},
        {"name": "p", "units": "Pa", "std_name": None, "provenance": None},
    ]


def test_normalize_applies_the_strong_form_provenance_unless_overridden():
    result = normalize_unknowns(
        ["u", {"name": "v"}, {"name": "p", "provenance": "human_edited"}],
        provenance="llm_derived",
    )
    assert [u["provenance"] for u in result] == [
        "llm_derived", "llm_derived", "human_edited"]


def test_normalize_rejects_bad_type():
    with pytest.raises(TypeError):
        normalize_unknowns([42])
