"""``lagged_coefficients`` (C42 round 5): which coefficients a sub-problem freezes, in time or in the Newton iteration."""

import copy

import pytest

pytest.importorskip("strictyaml")
pytest.importorskip("sympy")

from _formulation import TWO_PHASE, roundtrip, subproblem  # noqa: E402
from ymf.schema import load_ymf  # noqa: E402
from ymf.symbolic import split_problem  # noqa: E402
from ymf.vocabulary import LAG_WHEN, VocabularyError, check  # noqa: E402

DAMBREAK = TWO_PHASE / "dambreak.ymf"


def _lag(doc, lagged, sub="flow"):
    subproblem(doc, "level_set_five", sub)["lagged_coefficients"] = lagged
    return doc


def test_the_two_lags_are_different_words():
    assert set(LAG_WHEN) == {"previous_step", "previous_iterate"}


def test_a_lagged_state_coefficient_says_so_in_the_split_problem():
    doc = _lag(load_ymf(DAMBREAK), {"μ": "previous_iterate", "ρ": "previous_step"})
    flow = {s["name"]: s for s in split_problem(doc, "level_set_five")["subproblems"]}["flow"]
    assert flow["lagged_coefficients"] == {"μ": "previous_iterate", "ρ": "previous_step"}
    assert flow["coefficients"]["μ"]["lagged"] == "previous_iterate"
    assert flow["coefficients"]["ρ"]["lagged"] == "previous_step"
    assert "lagged" not in flow["coefficients"]["χ"]            # not stated: not lagged here


def test_what_is_not_stated_is_not_lagged():
    flow = {s["name"]: s for s in split_problem(load_ymf(DAMBREAK), "level_set_five")["subproblems"]}["flow"]
    assert "lagged_coefficients" not in flow
    assert all("lagged" not in c for c in flow["coefficients"].values())


@pytest.mark.parametrize("lagged, message", [
    ({"μ": "last_week"}, r"lagged_coefficients\.μ: lag when last_week: ymf knows previous_step, previous_iterate"),
    ({"viscosity": "previous_step"}, "viscosity is not a coefficient of the model; its coefficients are"),
    ({"φ": "previous_step"}, "φ is a field this sub-problem reads through lags, which already says when"),
])
def test_a_lagging_that_does_not_make_sense_is_refused_on_load(lagged, message):
    with pytest.raises(VocabularyError, match=message):
        check(_lag(load_ymf(DAMBREAK), lagged))


def test_lagging_survives_the_archive_and_is_part_of_the_key(tmp_path):
    from ymf import closure
    doc = _lag(load_ymf(DAMBREAK), {"μ": "previous_iterate"})
    path = tmp_path / "lagged.archive.ymf"
    closure.write(path, doc, {})
    again, _, _ = closure.load(path)
    assert subproblem(again, "level_set_five", "flow")["lagged_coefficients"] == {"μ": "previous_iterate"}
    realization = {"h": 0.0322}
    assert closure.input_digest(again, "level_set_five", realization) == closure.input_digest(
        doc, "level_set_five", realization)
    other = _lag(copy.deepcopy(doc), {"μ": "previous_step"})
    assert closure.input_digest(other, "level_set_five", realization) != closure.input_digest(
        doc, "level_set_five", realization)
