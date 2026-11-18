"""ymf 0.3.1: a spec written against ymf 0.2 (or main before 2026-11-17) loads unchanged.

The fixtures are proteus's own ``test/dispatch/dambreak.ymf`` and ``two_phase_level_set.ymf`` from ``cekees/main``
b9d98f5d, untouched (tests/fixtures/proteus_dispatch_0_2/README.md). ymf 0.3.0 refused them: RDLS's block states
``shock_capturing`` and no ``method``. They must load, mean what they said, and survive an archive.
"""

import pytest

pytest.importorskip("strictyaml")
pytest.importorskip("sympy")

from _formulation import FIXTURES, roundtrip, subproblem  # noqa: E402
from ymf import closure  # noqa: E402
from ymf.schema import load_ymf  # noqa: E402
from ymf.stabilization import missing  # noqa: E402
from ymf.symbolic import split_problem  # noqa: E402
from ymf.vocabulary import check  # noqa: E402

OLD = FIXTURES / "proteus_dispatch_0_2"
DAMBREAK = OLD / "dambreak.ymf"
SHOCK_ONLY = {"method": "none", "shock_capturing": {"factor": 0.9, "lagged": False}}


@pytest.fixture(scope="module")
def doc():
    return load_ymf(DAMBREAK)


def test_the_old_spec_loads_and_passes_the_checks(doc):
    check(doc)
    names = [e["name"] for e in doc["solution_paths"]["discretizations"]]
    assert names == ["level_set_five", "level_set_clsvof", "clsvof_projection"]


def test_the_old_blocks_mean_what_they_said(doc):
    split = {s["name"]: s for s in split_problem(doc, "level_set_five")["subproblems"]}
    # written `{method: ASGS, lag: true, shock_capturing: {factor: 0.5, lag: true}}`
    assert split["flow"]["stabilization"] == {"method": "ASGS", "lagged": True,
                                               "shock_capturing": {"factor": 0.5, "lagged": True}}
    # written `{method: SUPG, shock_capturing: {factor: 0.5, lag: true, sc_uref: 1.0, sc_beta: 1.5}}`
    for name in ("vof", "ncls"):
        assert split[name]["stabilization"] == {
            "method": "SUPG", "shock_capturing": {"factor": 0.5, "beta": 1.5, "reference": 1.0, "lagged": True}}
    # written `{shock_capturing: {factor: 0.9, lag: false}}`: shock capturing and no method; nothing was added
    assert split["rdls"]["stabilization"] == SHOCK_ONLY
    assert "stabilization" not in split["mcorr"]
    # the other two discretizations state none
    for name in ("level_set_clsvof", "clsvof_projection"):
        assert all("stabilization" not in s for s in split_problem(doc, name)["subproblems"])


def test_what_the_old_spec_leaves_unstated_is_still_reported_not_filled_in(doc):
    split = {s["name"]: s for s in split_problem(doc, "level_set_five")["subproblems"]}
    assert missing(split["rdls"]["stabilization"]) == ["shock_capturing.beta", "shock_capturing.reference"]
    assert missing(split["vof"]["stabilization"]) == ["tau", "lagged"]
    assert missing(split["flow"]["stabilization"]) == ["tau", "shock_capturing.beta", "shock_capturing.reference"]


def test_the_old_spec_round_trips_through_an_archive(tmp_path):
    realization = {"h": 0.0322}
    spec, again, digest = roundtrip(tmp_path, DAMBREAK, "level_set_five", realization)
    # the key is the old spec's: reading it again from its archive gives the same one, and the archived blocks resolve
    # to what the fixture's did
    for name in ("flow", "vof", "ncls", "rdls", "mcorr"):
        assert subproblem(again, "level_set_five", name).get("stabilization") == subproblem(
            spec, "level_set_five", name).get("stabilization")
    split = {s["name"]: s for s in split_problem(again, "level_set_five")["subproblems"]}
    assert split["rdls"]["stabilization"] == SHOCK_ONLY
    # and the statement is part of the key: the same spec with RDLS's shock capturing changed has another one
    other = load_ymf(DAMBREAK)
    subproblem(other, "level_set_five", "rdls")["stabilization"] = {"shock_capturing": {"factor": 0.8, "lag": False}}
    assert closure.input_digest(other, "level_set_five", realization) != digest


def test_the_old_block_and_the_same_in_ymfs_new_words_resolve_alike(doc):
    new = load_ymf(DAMBREAK)
    subproblem(new, "level_set_five", "rdls")["stabilization"] = {
        "method": "none", "shock_capturing": {"factor": 0.9, "lagged": False}}
    rdls = lambda d: {s["name"]: s for s in split_problem(d, "level_set_five")["subproblems"]}["rdls"]  # noqa: E731
    assert rdls(doc)["stabilization"] == rdls(new)["stabilization"] == SHOCK_ONLY
