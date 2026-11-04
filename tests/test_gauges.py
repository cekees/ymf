"""Gauges: what a run measures as it goes.

The shape FORCE_HISTORY_GAUGES_PLAN step 4 planned (C22): point, line and
line-integral gauges of named fields, and the boundary-integral gauge of
forces and moments per boundary region, which replaces RANS2P's and
RANS3PF's history files.
"""

import copy
from pathlib import Path

import pytest

pytest.importorskip("strictyaml")

from ymf.references import YmfReferenceError, resolve  # noqa: E402
from ymf.schema import load_ymf  # noqa: E402

TWO_PHASE = Path(__file__).parent / "fixtures" / "two_phase"


@pytest.fixture
def dambreak():
    return load_ymf(TWO_PHASE / "dambreak.ymf")


def _with(doc, gauge):
    doc = copy.deepcopy(doc)
    doc["model"]["gauges"] = [gauge]
    return doc


def test_the_dam_break_carries_its_decks_gauges(dambreak):
    gauges = {g["name"]: g for g in dambreak["model"]["gauges"]}
    assert len(gauges["heights"]["lines"]) == 4
    assert gauges["pressure"]["points"][2] == [3.22, 0.12]
    assert gauges["pressure"]["file"] == "pressure.csv"


def test_a_boundary_integral_gauge(dambreak):
    resolve(_with(dambreak, {"kind": "boundary_integral", "tags": ["walls"],
                             "quantities": ["force", "moment"],
                             "moment_point": [1.61, 0.9], "sample_rate": 0.01,
                             "file": "wall_forces.csv"}))


@pytest.mark.parametrize("gauge, message", [
    ({"kind": "point", "fields": ["p"]}, "a point gauge needs points"),
    ({"kind": "line", "points": [[0.0, 0.0]], "fields": ["p"], "lines": [[[0, 0], [1, 1]]]},
     "a line gauge takes no points"),
    ({"kind": "boundary_integral", "tags": ["walls"]}, "needs quantities"),
    ({"kind": "boundary_integral", "tags": ["walls"], "quantities": ["force"],
      "fields": ["p"]}, "takes no fields"),
    ({"kind": "point", "fields": ["q"], "points": [[0.0, 0.0]]},
     "names 'q', which is not an unknown or coefficient"),
    ({"kind": "boundary_integral", "tags": ["lid"], "quantities": ["force"]},
     "names 'lid', which is not a boundary region"),
    ({"kind": "line_integral", "fields": ["φ"], "lines": [[[0, 0], [1, 1], [2, 2]]]},
     "a line is its two ends, not 3 points"),
    ({"kind": "point", "fields": ["p"], "points": [[0.0, 0.0], [1.0, 1.0, 1.0]]},
     "do not all have the same dimension"),
])
def test_what_a_gauge_needs_is_checked(dambreak, gauge, message):
    with pytest.raises(YmfReferenceError, match=message):
        resolve(_with(dambreak, gauge))


def test_a_gauge_of_a_coefficient_of_the_state(dambreak):
    resolve(_with(dambreak, {"kind": "point", "fields": ["ρ"], "points": [[1.0, 0.5]]}))
