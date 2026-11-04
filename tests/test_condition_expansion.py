"""Boundary-condition aliases expand into the abstract forms (cekees, 2026-11-04).

slip, open, inflow and characteristic, and a model's condition sets, are
shorthand; what is solved, and what a closed spec records, is Dirichlet,
Neumann and Robin conditions on a variable or a projection of it, with their
values and with what they expanded from. ymf.conditions holds the meaning as
data (ALIASES) and expands a document (expand).
"""

import copy
from pathlib import Path

import pytest

from ymf.conditions import ALIASES, expand, expand_condition

FIXTURES = Path(__file__).parent / "fixtures"


@pytest.mark.parametrize("alias, rank, given, expected", [
    ("slip", 1, {}, [
        {"type": "dirichlet", "component": "normal", "value": 0.0},
        {"type": "neumann", "component": "tangential", "value": 0.0}]),
    ("slip", 1, {"formula": "(t, 0)"}, [
        {"type": "dirichlet", "component": "normal", "formula": "(t, 0)"},
        {"type": "neumann", "component": "tangential", "value": 0.0}]),
    ("open", 0, {}, [{"type": "neumann", "value": 0.0}]),
    ("inflow", 0, {"value": 1.0, "units": "1"}, [
        {"type": "dirichlet", "only": "inflow", "value": 1.0, "units": "1"}]),
    ("characteristic", 0, {"series": [[0.0, 0.1], [5.0, 0.2]]}, [
        {"type": "dirichlet", "only": "incoming_characteristics",
         "series": [[0.0, 0.1], [5.0, 0.2]]}]),
])
def test_each_alias_expands_to_its_documented_components(alias, rank, given, expected):
    bc = dict({"region": "r", "variable": "u", "type": alias}, **given)
    parts = expand_condition(bc, rank)
    assert [{k: v for k, v in p.items() if k not in ("region", "variable", "expanded_from")}
            for p in parts] == expected
    assert all(p["region"] == "r" and p["variable"] == "u" for p in parts)
    assert all(p["expanded_from"] == {"type": alias} for p in parts)


def test_every_alias_is_made_of_abstract_forms_only():
    for name, alias in ALIASES.items():
        assert alias["means"]
        assert {c["type"] for c in alias["components"]} <= {"dirichlet", "neumann", "robin"}


def test_an_abstract_condition_is_itself():
    bc = {"region": "r", "variable": "u", "type": "robin", "alpha": "1", "beta": "2",
          "value": 0.5}
    assert expand_condition(bc) == [bc]


def test_slip_needs_a_vector():
    with pytest.raises(ValueError, match="a slip condition is on a vector, and h is not one"):
        expand_condition({"region": "r", "variable": "h", "type": "slip"}, 0)


def _strong(**kw):
    return {"model": {"name": "m", "strong_form": dict(
        {"unknowns": [{"name": "v", "rank": 1}, "p"], "domain": "Ω = [0, 1] × [0, 1]",
         "boundary_regions": [{"name": "walls", "where": "y = 0"},
                              {"name": "top", "where": "y = 1"}]}, **kw)}}


def test_a_document_without_shorthand_comes_back_unchanged():
    doc = _strong(boundary_conditions=[{"region": "walls", "variable": "v",
                                        "type": "dirichlet", "value": 0.0}])
    assert expand(doc) is doc


def test_sets_and_aliases_expand_in_order_and_the_sets_go():
    doc = _strong(
        boundary_conditions=[{"region": "walls", "variable": "p", "type": "neumann",
                              "value": 0.0}],
        condition_sets={"atmosphere": [{"variable": "p", "type": "dirichlet", "value": 0.0},
                                       {"variable": "v", "type": "open"}]},
        boundary_sets=[{"region": "top", "set": "atmosphere"}])
    strong = expand(doc)["model"]["strong_form"]
    assert "condition_sets" not in strong and "boundary_sets" not in strong
    assert strong["boundary_conditions"] == [
        {"region": "walls", "variable": "p", "type": "neumann", "value": 0.0},
        {"region": "top", "variable": "p", "type": "dirichlet", "value": 0.0,
         "expanded_from": {"set": "atmosphere"}},
        {"region": "top", "variable": "v", "type": "neumann", "value": 0.0,
         "expanded_from": {"type": "open", "set": "atmosphere"}}]


def test_expanding_is_idempotent_and_leaves_its_input_alone():
    doc = _strong(boundary_conditions=[{"region": "walls", "variable": "v", "type": "slip"}])
    before = copy.deepcopy(doc)
    once = expand(doc)
    assert doc == before
    assert expand(once) is once


def test_what_does_not_expand_is_refused():
    with pytest.raises(ValueError, match="no condition set named 'nowhere'"):
        expand(_strong(boundary_sets=[{"region": "top", "set": "nowhere"}]))
