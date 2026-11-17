"""Shared by the round-5 formulation tests (stabilization, lagging, boundary kinds, dispersion, cells): the fixtures,
and a round trip through an archive."""

from pathlib import Path

FIXTURES = Path(__file__).parent / "fixtures"
TWO_PHASE = FIXTURES / "two_phase"
SHALLOW_WATER = FIXTURES / "shallow_water"


def roundtrip(tmp_path, spec_path, discretization, realization):
    """Load a spec, write its archive, read that back, and once more: ``(spec, again, digest)`` where the digest (the
    input's key) is the same each time, so that a word a spec states is part of what an archive records and not lost
    or rewritten by it."""
    from ymf import closure
    spec, _, _ = closure.load(spec_path)
    digest = closure.input_digest(spec, discretization, realization)
    archive = tmp_path / "round.archive.ymf"
    closure.write(archive, spec, {})
    again, _, _ = closure.load(archive)
    assert closure.input_digest(again, discretization, realization) == digest
    closure.write(archive, again, {})
    third, _, _ = closure.load(archive)
    assert closure.input_digest(third, discretization, realization) == digest
    return spec, again, digest


def subproblem(spec, discretization, name):
    entry = next(e for e in spec["solution_paths"]["discretizations"] if e["name"] == discretization)
    return next(s for s in entry["subproblems"] if s["name"] == name)
