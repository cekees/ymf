"""A boundary value taken from another sub-problem's trace (C42 round 7; C19 gap note G8).

NCLS's inflow penalty compares its level set with RDLS's face values: the condition's value is not an expression but the
boundary values of another sub-problem's field. ``trace: {of, via, at}`` says it, in place of value, formula or series.
"""

import copy

import pytest

pytest.importorskip("strictyaml")
pytest.importorskip("sympy")

from _formulation import TWO_PHASE, roundtrip  # noqa: E402
from ymf import closure  # noqa: E402
from ymf.references import YmfReferenceError, resolve  # noqa: E402
from ymf.symbolic import SymbolicError, split_problem  # noqa: E402

FORMULATION = TWO_PHASE / "dambreak_formulation.ymf"
REALIZATION = {"h": 0.0322}
TRACE = {"region": "top", "variable": "φ", "type": "dirichlet", "imposed": "weak",
         "trace": {"of": "φ", "via": "rdls", "at": "previous_step"}}


def _spec(*conditions):
    spec, _, _ = closure.load(FORMULATION)
    strong = spec["model"]["strong_form"]
    strong["boundary_conditions"] = [copy.deepcopy(c) for c in conditions]
    return spec


def _sub(spec, name):
    return next(s for s in split_problem(spec, "level_set_five")["subproblems"] if s["name"] == name)


def test_a_trace_reaches_the_sub_problem_that_solves_the_field():
    spec = _spec(TRACE)
    resolve(spec)
    (entry,) = [c for c in _sub(spec, "ncls")["boundary_conditions"] if c["variable"] == "φ"]
    assert entry["trace"] == {"of": "φ", "via": "rdls", "at": "previous_step"}
    assert "value" not in entry and "series" not in entry and entry["kind"] == "weak"
    # the other sub-problems do not impose it
    assert not [c for name in ("flow", "vof", "rdls") for c in _sub(spec, name)["boundary_conditions"]
                if c.get("trace")]
    assert not [c for c in _sub(spec, "rdls")["boundary_conditions"] if c["variable"] == "φ"]


def test_of_defaults_to_the_conditions_variable_and_a_current_step_read_needs_an_earlier_sub_problem():
    cond = dict(TRACE, trace={"via": "flow", "at": "current_step"})      # flow is solved before ncls
    (entry,) = [c for c in _sub(_spec(cond), "ncls")["boundary_conditions"] if c["variable"] == "φ"]
    assert entry["trace"] == {"of": "φ", "via": "flow", "at": "current_step"}
    late = dict(TRACE, trace={"via": "rdls", "at": "current_step"})      # rdls is solved after ncls
    with pytest.raises(SymbolicError, match="solved after ncls in the step; say previous_step"):
        _sub(_spec(late), "ncls")


def test_the_sub_problem_whose_trace_it_is_imposes_nothing():
    # rdls also solves φ: the condition is for the sub-problems that read rdls's face values, not for rdls
    assert not [c for c in _sub(_spec(TRACE), "rdls")["boundary_conditions"] if c["variable"] == "φ"]


@pytest.mark.parametrize("condition, message", [
    (dict(TRACE, value=1.0), "takes no value, formula or series"),
    (dict(TRACE, formula="y"), "takes no value, formula or series"),
    (dict(TRACE, type="neumann"), "a neumann condition takes no trace"),
    (dict(TRACE, trace={"via": "nowhere", "at": "previous_step"}), "'nowhere'.*sub-problem"),
    (dict(TRACE, trace={"via": "rdls", "of": "ψ", "at": "previous_step"}), "'ψ'.*model unknown"),
    ({"region": "top", "variable": "v", "type": "dirichlet", "component": "normal",
      "trace": {"via": "rdls", "at": "previous_step"}}, "takes no component"),
])
def test_what_a_trace_needs_is_checked(condition, message):
    with pytest.raises(YmfReferenceError, match=message):
        resolve(_spec(condition))


def test_one_adr_system_runs_no_trace():
    from ymf.schema import load_ymf
    from pathlib import Path
    from ymf.symbolic import adr_problem
    doc = load_ymf(Path(__file__).parent.parent / "examples" / "poisson.ymf")
    bc = doc["model"]["strong_form"]["boundary_conditions"][0]
    bc.pop("value", None)
    bc.pop("formula", None)
    bc["trace"] = {"via": "x", "at": "previous_step"}
    with pytest.raises(SymbolicError, match="another sub-problem's trace as its value runs in a split"):
        adr_problem(doc)


def test_a_trace_round_trips_through_an_archive_and_is_part_of_the_key(tmp_path, monkeypatch):
    spec = _spec(TRACE)
    path = tmp_path / "traced.ymf"
    closure.write(path, spec, {})
    # an extending file would be the usual form; the closed spec is enough to show the word survives and keys
    again, _, _ = closure.load(path)
    assert again["model"]["strong_form"]["boundary_conditions"][0]["trace"] == TRACE["trace"]
    assert _sub(again, "ncls")["boundary_conditions"] == _sub(spec, "ncls")["boundary_conditions"]
    digest = closure.input_digest(spec, "level_set_five", REALIZATION)
    assert closure.input_digest(again, "level_set_five", REALIZATION) == digest
    other = _spec(dict(TRACE, trace={"via": "flow", "at": "current_step"}))
    assert closure.input_digest(other, "level_set_five", REALIZATION) != digest
