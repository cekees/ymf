"""``extends`` pairs initial conditions by field (C5 finding Y3): a file
that restates one field's initial condition replaces that one, and inherits
the rest."""

from pathlib import Path

import pytest

pytest.importorskip("strictyaml")

from ymf.schema import load_ymf  # noqa: E402

TWO_PHASE = Path(__file__).parent / "fixtures" / "two_phase"


def _ics(doc):
    return {ic["field"]: ic for ic in doc["model"]["strong_form"]["initial_conditions"]}


def test_one_field_restated_the_others_inherited(tmp_path):
    for name in ("two_phase_level_set.ymf", "dambreak.ymf"):
        (tmp_path / name).write_text((TWO_PHASE / name).read_text(encoding="utf-8"),
                                     encoding="utf-8")
    (tmp_path / "column.ymf").write_text(
        "extends: dambreak.ymf\n"
        "model:\n"
        "  strong_form:\n"
        "    initial_conditions:\n"
        "      - {field: φ, type: function, formula: \"y − 0.6\", units: m}\n"
        "      - {field: p, type: function, formula: \"998.2*9.8*(0.6 − y)\"}\n",
        encoding="utf-8")
    doc = load_ymf(tmp_path / "column.ymf")
    ics = _ics(doc)
    assert list(ics) == ["v", "p", "φ"]                          # the parent's order
    assert ics["v"] == {"field": "v", "type": "constant", "value": 0.0}
    # replaced whole: a constant's value does not survive into a function
    assert ics["p"] == {"field": "p", "type": "function", "formula": "998.2*9.8*(0.6 − y)"}
    assert ics["φ"]["formula"] == "y − 0.6"
    overrides = doc["composition"]["overrides"]
    assert {o["path"] for o in overrides} >= {"model.strong_form.initial_conditions[φ]",
                                              "model.strong_form.initial_conditions[p]"}


def test_a_new_field_is_appended_and_one_can_be_removed(tmp_path):
    for name in ("two_phase_level_set.ymf", "dambreak.ymf"):
        (tmp_path / name).write_text((TWO_PHASE / name).read_text(encoding="utf-8"),
                                     encoding="utf-8")
    (tmp_path / "drop.ymf").write_text(
        "extends: dambreak.ymf\n"
        "ymf_state: abstract\n"
        "remove:\n"
        "  - \"model.strong_form.initial_conditions[p]\"\n",
        encoding="utf-8")
    doc = load_ymf(tmp_path / "drop.ymf")
    assert list(_ics(doc)) == ["v", "φ"]
