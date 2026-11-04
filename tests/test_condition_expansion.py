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


# --- in closure -------------------------------------------------------------------------

strictyaml = pytest.importorskip("strictyaml")
from ymf import closure  # noqa: E402
from ymf.schema import load_ymf  # noqa: E402


def _keys(doc):
    return sorted(k for k, _, _ in closure.planned_outputs(doc))


@pytest.mark.parametrize("fixture", ["two_phase/dambreak.ymf", "shallow_water/torino.ymf",
                                     "shallow_water/lake_at_rest.ymf"])
def test_a_closed_spec_holds_the_expansion_and_round_trips(fixture, tmp_path):
    spec = load_ymf(FIXTURES / fixture)
    archive = tmp_path / "closed.archive.ymf"
    closure.write(archive, spec, {})
    reloaded, outputs, dropped = closure.load(archive)
    strong = reloaded["model"]["strong_form"]
    assert "boundary_sets" not in strong and "condition_sets" not in strong
    assert all(bc["type"] in ("dirichlet", "neumann", "robin", "periodic")
               for bc in strong["boundary_conditions"])
    assert strong["boundary_conditions"] == expand(spec)["model"]["strong_form"][
        "boundary_conditions"]
    assert _keys(reloaded) == _keys(spec)          # the spec and its closure key alike
    closure.write(tmp_path / "again.archive.ymf", reloaded, {})
    assert (tmp_path / "again.archive.ymf").read_text() == archive.read_text()


def _without_provenance(doc):
    doc = copy.deepcopy(doc)
    for bc in doc["model"]["strong_form"]["boundary_conditions"]:
        bc.pop("expanded_from", None)
    return doc


def test_aliases_and_the_same_conditions_written_out_close_alike(tmp_path):
    """The dam break with its sets, and with what they mean written by hand."""
    aliased = load_ymf(FIXTURES / "two_phase" / "dambreak.ymf")
    (tmp_path / "two_phase_level_set.ymf").write_text(
        (FIXTURES / "two_phase" / "two_phase_level_set.ymf").read_text(encoding="utf-8")
        .split("    # what the cases' walls and tops mean")[0], encoding="utf-8")
    text = (FIXTURES / "two_phase" / "dambreak.ymf").read_text(encoding="utf-8")
    text = text.replace('''    boundary_sets:
      - {region: walls, set: free_slip}
      - {region: top, set: atmosphere}''', '''    boundary_conditions:
      - {region: walls, variable: v, type: dirichlet, component: normal, value: 0.0}
      - {region: walls, variable: v, type: neumann, component: tangential, value: 0.0}
      - {region: top, variable: p, type: dirichlet, value: 0.0}
      - {region: top, variable: v, type: neumann, value: 0.0}
      - {region: top, variable: χ, type: dirichlet, only: inflow, value: 1.0}''')
    (tmp_path / "dambreak.ymf").write_text(text, encoding="utf-8")
    by_hand = load_ymf(tmp_path / "dambreak.ymf")
    a = _without_provenance(expand(aliased))
    b = _without_provenance(expand(by_hand))
    a.pop("composition", None)
    b.pop("composition", None)
    assert a == b
    closure.write(tmp_path / "a.archive.ymf", aliased, {})
    closure.write(tmp_path / "b.archive.ymf", by_hand, {})
    closed = []
    for name in ("a", "b"):
        doc, _, _ = closure.load(tmp_path / ("%s.archive.ymf" % name))
        doc = _without_provenance(doc)
        doc.pop("composition", None)
        closed.append(doc)
    assert closed[0] == closed[1]            # the two closed documents, provenance aside


def test_a_projection_is_of_a_vector():
    from ymf.references import YmfReferenceError, resolve
    doc = load_ymf(FIXTURES.parent.parent / "examples" / "heat_equation.ymf")
    doc["model"]["strong_form"]["boundary_conditions"][0].update(component="normal")
    with pytest.raises(YmfReferenceError, match="a normal component is a vector's, and T is "
                                                "a scalar"):
        resolve(doc)


def test_a_reader_sees_only_the_expansion():
    """split_problem and adr_problem expand first; units checks the expansion."""
    pytest.importorskip("sympy")
    from ymf.symbolic import SymbolicError, adr_problem
    doc = load_ymf(FIXTURES / "shallow_water" / "lake_at_rest.ymf")
    with pytest.raises(SymbolicError, match="normal component \\(expanded from slip\\)"):
        adr_problem(doc)
