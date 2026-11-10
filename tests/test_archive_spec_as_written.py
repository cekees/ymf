"""An archive's specification part reads back as written (C5 finding Y1).

The archive is read with pyyaml, which typed a free-form value: ``1.0e6``
under a ``proteus:`` hint came back as the float 1e6, was written again as
``1000000.0``, and strictyaml read that string -- another input than the
spec's ``'1.0e6'``, so the archive's outputs no longer matched their key and
``--check`` found nothing to check."""

from pathlib import Path

import pytest

pytest.importorskip("strictyaml")
pytest.importorskip("yaml")

from ymf import closure  # noqa: E402

TWO_PHASE = Path(__file__).parent / "fixtures" / "two_phase"

HINTS = "proteus: {model: RANS2P, useVF: 1.0, eb_penalty_constant: 1.0e6, tol: 1e-3, " \
        "on: true, mask: '0x10', n: 007}"


def _write(tmp_path):
    for name in ("two_phase_level_set.ymf", "dambreak.ymf"):
        text = (TWO_PHASE / name).read_text(encoding="utf-8")
        text = text.replace("proteus: {model: RANS2P, useVF: 1.0}", HINTS)
        (tmp_path / name).write_text(text, encoding="utf-8")
    return tmp_path / "dambreak.ymf"


def _hints(spec):
    entry = spec["solution_paths"]["discretizations"][0]
    return entry["subproblems"][0]["proteus"]


def test_free_form_values_and_the_key_survive_the_archive(tmp_path):
    spec, _, _ = closure.load(_write(tmp_path))
    realization = {"h": 0.0322}
    before = closure.input_digest(spec, "level_set_five", realization)
    assert _hints(spec)["eb_penalty_constant"] == "1.0e6"
    archive = tmp_path / "dambreak.archive.ymf"
    closure.write(archive, spec, {})
    again, _, _ = closure.load(archive)
    assert _hints(again) == _hints(spec)
    assert closure.input_digest(again, "level_set_five", realization) == before
    # and once more: an archive of the archive is the same input
    closure.write(archive, again, {})
    third, _, _ = closure.load(archive)
    assert closure.input_digest(third, "level_set_five", realization) == before


def test_typed_values_keep_their_types(tmp_path):
    spec, _, _ = closure.load(_write(tmp_path))
    archive = tmp_path / "dambreak.archive.ymf"
    closure.write(archive, spec, {})
    again, _, _ = closure.load(archive)
    coefficients = again["model"]["strong_form"]["coefficients"]
    assert coefficients["νw"]["value"] == 1.004e-6
    entry = again["solution_paths"]["discretizations"][0]
    assert entry["mesh"]["h"] == 0.0322 and entry["step_control"]["step_exact"] is True
    assert entry["subproblems"][0]["solver"]["max_iterations"] == 100
