"""ymf.stabilization (C42 round 5): the stabilization and shock-capturing words, their parameters, checked on load and
carried through the split problem and the archive."""

import pytest

pytest.importorskip("strictyaml")
sympy = pytest.importorskip("sympy")

from _formulation import TWO_PHASE, roundtrip, subproblem  # noqa: E402
from ymf import stabilization  # noqa: E402
from ymf.references import YmfReferenceError  # noqa: E402
from ymf.schema import load_ymf  # noqa: E402
from ymf.stabilization import StabilizationError, missing, resolve  # noqa: E402
from ymf.symbolic import split_problem  # noqa: E402
from ymf.vocabulary import VocabularyError, canonical, check  # noqa: E402

FORMULATION = TWO_PHASE / "dambreak_formulation.ymf"
DAMBREAK = TWO_PHASE / "dambreak.ymf"


# --- resolve -----------------------------------------------------------------------------------------------------

def test_a_bare_word_is_a_block_and_nothing_is_added():
    assert resolve("SUPG") == {"method": "SUPG"}
    assert resolve({"method": "SUPG", "tau": "element_size"}) == {"method": "SUPG", "tau": "element_size"}


def test_the_keys_the_dispatcher_reads_are_aliases_of_ymfs():
    old = {"method": "SUPG", "lag": True, "shock_capturing": {"factor": 0.5, "sc_beta": 1.5, "sc_uref": 1.0, "lag": True}}
    new = {"method": "SUPG", "lagged": True,
           "shock_capturing": {"factor": 0.5, "beta": 1.5, "reference": 1.0, "lagged": True}}
    assert resolve(old) == resolve(new) == new
    assert canonical("shock_capturing_key", "sc_uref") == "reference"
    assert resolve({"method": "shock_capturing"}) == {"method": "residual_shock_capturing"}


def test_shock_capturing_none_is_stated_and_is_not_absent():
    assert resolve({"method": "SUPG", "shock_capturing": "none"}) == {"method": "SUPG", "shock_capturing": None}
    assert "shock_capturing" not in resolve("SUPG")


@pytest.mark.parametrize("block, message", [
    ({"tau": "metric"}, "say the method"),
    ({}, "say the method"),
    ({"tau": "metric", "shock_capturing": {"factor": 0.9}},
     "no method is stated.*tau is read by a method: say which"),
    ({"tau": "metric", "lagged": True, "shock_capturing": {"factor": 0.9}},
     "lagged, tau are read by a method"),
    ({"method": "none", "tau": "metric"}, "none does not read tau; its keys are method, shock_capturing"),
    ({"method": "none", "lagged": True}, "none does not read lagged"),
    ({"shock_capturing": {"beta": 1.5}}, "factor is required"),
    ({"method": "GLS"}, "stabilization GLS: ymf knows"),
    ({"method": "SUPG", "tau": "h"}, "tau: h is not one of metric, element_size"),
    ({"method": "SUPG", "lagged": "maybe"}, "'maybe' is not true or false"),
    ({"method": "SUPG", "lagged": True, "lag": False}, "lagged is stated twice"),
    ({"method": "SUPG", "shock_capturing": {"factor": 0.5, "sc_beta": 1.2, "beta": 1.2}}, "beta is stated twice"),
    ({"method": "SUPG", "shock_capturing": {"beta": 1.5}}, "factor is required"),
    ({"method": "SUPG", "shock_capturing": {"factor": -1.0}}, "factor: -1 is negative"),
    ({"method": "SUPG", "shock_capturing": {"factor": 0.5, "beta": 0.5}}, "beta: 0.5 is outside \\[1, 2\\]"),
    ({"method": "SUPG", "shock_capturing": {"factor": 0.5, "beta": 2.5}}, "beta: 2.5 is outside \\[1, 2\\]"),
    ({"method": "SUPG", "shock_capturing": {"factor": 0.5, "reference": 0.0}}, "reference: 0 is not positive"),
    ({"method": "SUPG", "shock_capturing": {"factor": 0.5, "gamma": 1.0}}, "no key gamma in shock_capturing"),
    ({"method": "SUPG", "shock_capturing": "some"}, "'some' is not a block"),
    ({"method": "entropy_viscosity", "tau": "metric"},
     "entropy_viscosity does not read tau; it has no keys here .*scheme's options"),
    ({"method": "graph_viscosity_smoothness", "shock_capturing": {"factor": 1.0}},
     "does not read shock_capturing"),
    ({"method": "residual_shock_capturing", "lagged": True},
     "residual_shock_capturing does not read lagged; its keys are method, tau, shock_capturing"),
])
def test_a_stabilization_that_does_not_make_sense_is_refused(block, message):
    with pytest.raises(StabilizationError, match=message):
        resolve(block)


# --- no method: none, and shock capturing alone (0.3.1) ------------------------------------------------------------

def test_none_says_there_is_no_stabilization_and_adds_nothing():
    assert resolve("none") == {"method": "none"} == resolve({"method": "none"})
    assert canonical("stabilization", "none") == "none"
    # none with shock capturing stated none is the same statement, spelled out
    assert resolve({"method": "none", "shock_capturing": "none"}) == {"method": "none", "shock_capturing": None}


def test_a_block_that_states_only_shock_capturing_is_shock_capturing_with_no_method():
    # proteus's RDLS block in test/dispatch at cekees/main b9d98f5d (ymf 0.2: a free map; 0.3.0: refused for want of method)
    old = {"shock_capturing": {"factor": 0.9, "lag": False}}
    assert resolve(old) == {"method": "none", "shock_capturing": {"factor": 0.9, "lagged": False}}
    assert resolve(old) == resolve({"method": "none", "shock_capturing": {"factor": 0.9, "lagged": False}})
    # not a method ymf picks: it is not residual_shock_capturing (a method that also takes tau), nor SUPG
    assert resolve(old)["method"] not in ("SUPG", "residual_shock_capturing")
    assert resolve({"shock_capturing": "none"}) == {"method": "none", "shock_capturing": None}


def test_missing_asks_none_only_for_a_shock_capturing_block_it_states():
    assert missing("none") == [] and missing({"method": "none", "shock_capturing": "none"}) == []
    assert missing({"shock_capturing": {"factor": 0.9, "lag": False}}) == [
        "shock_capturing.beta", "shock_capturing.reference"]
    assert missing({"shock_capturing": {"factor": 0.9, "beta": 1.0, "reference": 1.0, "lagged": False}}) == []
    # nothing was filled in on the way: the same two keys are still unstated
    assert "beta" not in resolve({"shock_capturing": {"factor": 0.9}})["shock_capturing"]
    # SUPG still asks for everything it reads, shock capturing stated or not
    assert missing("SUPG") == ["tau", "lagged", "shock_capturing"]


def test_a_spec_states_none_and_shock_capturing_alone_and_split_problem_carries_them():
    doc = load_ymf(FORMULATION)
    _state(doc, "mcorr", "none")
    _state(doc, "rdls", {"shock_capturing": {"factor": 0.9, "lag": False}})
    check(doc)
    split = {s["name"]: s for s in split_problem(doc, "level_set_five")["subproblems"]}
    assert split["mcorr"]["stabilization"] == {"method": "none"}
    assert split["rdls"]["stabilization"] == {"method": "none", "shock_capturing": {"factor": 0.9, "lagged": False}}
    # an entry's none is the default of the sub-problems that state nothing; a sub-problem that states one replaces it
    entry = doc["solution_paths"]["discretizations"][0]
    entry["stabilization"] = {"method": "none"}
    del subproblem(doc, "level_set_five", "vof")["stabilization"]
    split = {s["name"]: s for s in split_problem(doc, "level_set_five")["subproblems"]}
    assert split["vof"]["stabilization"] == {"method": "none"}
    assert split["flow"]["stabilization"]["method"] == "ASGS"


def test_a_block_that_states_shock_capturing_with_a_key_only_a_method_reads_is_refused_on_load():
    doc = _state(load_ymf(FORMULATION), "rdls", {"tau": "metric", "shock_capturing": {"factor": 0.9}})
    with pytest.raises(VocabularyError, match=r"subproblems\[rdls\]\.stabilization: no method is stated"):
        check(doc)
    doc = _state(load_ymf(FORMULATION), "rdls", {"tau": "metric"})
    with pytest.raises(VocabularyError, match="say the method"):
        check(doc)


# --- what a generator still needs stated ------------------------------------------------------------------------

def test_missing_says_what_the_generator_cannot_guess():
    assert missing("SUPG") == ["tau", "lagged", "shock_capturing"]
    assert missing({"method": "SUPG", "tau": "metric", "lagged": False,
                    "shock_capturing": {"factor": 0.5, "lagged": True}}) == [
        "shock_capturing.beta", "shock_capturing.reference"]
    assert missing({"method": "SUPG", "tau": "metric", "lagged": False, "shock_capturing": "none"}) == []
    assert missing("entropy_viscosity") == []          # its parameters are the scheme's


def test_reference_is_not_missing_at_beta_1_and_still_valid_when_stated():
    rdls = {"method": "SUPG", "tau": "metric", "lagged": False,
            "shock_capturing": {"factor": 0.9, "beta": 1.0, "lagged": False}}      # proteus's RDLS
    assert missing(rdls) == []
    assert missing({"shock_capturing": {"factor": 0.9, "beta": 1, "lag": False}}) == []      # also the 0.2 spelling
    assert missing({"shock_capturing": {"factor": 0.9, "sc_beta": 1.0, "lag": False}}) == []
    # beta other than 1, or not stated, still needs it
    assert missing({"shock_capturing": {"factor": 0.9, "beta": 1.5, "lagged": False}}) == ["shock_capturing.reference"]
    assert missing({"shock_capturing": {"factor": 0.9, "lagged": False}}) == [
        "shock_capturing.beta", "shock_capturing.reference"]
    # stating it is still fine, and still checked: positive
    stated = dict(rdls, shock_capturing=dict(rdls["shock_capturing"], reference=2.0))
    assert missing(stated) == [] and resolve(stated)["shock_capturing"]["reference"] == 2.0
    with pytest.raises(StabilizationError, match="reference"):
        resolve(dict(rdls, shock_capturing=dict(rdls["shock_capturing"], reference=-1.0)))


# --- the spec ----------------------------------------------------------------------------------------------------

@pytest.fixture(scope="module")
def formulation():
    return load_ymf(FORMULATION)


def test_the_dam_break_states_the_stabilization_the_dispatcher_used_to_assume(formulation):
    split = {s["name"]: s for s in split_problem(formulation, "level_set_five")["subproblems"]}
    # proteus.DispatchTable at cekees/main 92f24ad9: RANS2P ASGS, subgrid lag, shock capturing 0.5 lagged; VOF and NCLS
    # SUPG with 0.5, sc_beta 1.5, sc_uref 1.0, lagged; RDLS SUPG with 0.9 not lagged; MCorr none
    assert split["flow"]["stabilization"] == {"method": "ASGS", "tau": "metric", "lagged": True,
                                               "shock_capturing": {"factor": 0.5, "lagged": True}}
    for name in ("vof", "ncls"):
        assert split[name]["stabilization"] == {
            "method": "SUPG", "tau": "metric", "lagged": False,
            "shock_capturing": {"factor": 0.5, "beta": 1.5, "reference": 1.0, "lagged": True}}
    assert split["rdls"]["stabilization"]["shock_capturing"] == {"factor": 0.9, "lagged": False}
    assert "stabilization" not in split["mcorr"]


def test_vof_and_ncls_state_what_the_generator_spike_declared_by_hand(formulation):
    # generator/make_inputs.py (proteus#298): STAB = {"method": "SUPG", "tau": "metric", "shock_capturing": {"lagged": True}}
    declared = {"method": "SUPG", "tau": "metric", "shock_capturing": {"lagged": True}}
    split = {s["name"]: s for s in split_problem(formulation, "level_set_five")["subproblems"]}
    for name in ("vof", "ncls"):
        stated = split[name]["stabilization"]
        assert stated["method"] == declared["method"] and stated["tau"] == declared["tau"]
        assert stated["shock_capturing"]["lagged"] is declared["shock_capturing"]["lagged"]
        # and what it needed beyond that, which make_inputs.py could not say
        assert missing(stated) == []


def test_the_dam_break_without_the_formulation_states_none_of_it():
    split = {s["name"]: s for s in split_problem(load_ymf(DAMBREAK), "level_set_five")["subproblems"]}
    assert all("stabilization" not in s for s in split.values())


def test_an_entry_states_the_default_of_its_sub_problems_and_a_sub_problem_replaces_it_whole(tmp_path):
    doc = load_ymf(FORMULATION)
    entry = doc["solution_paths"]["discretizations"][0]
    entry["stabilization"] = {"method": "SUPG", "tau": "element_size"}
    for sub in entry["subproblems"]:
        if sub["name"] in ("vof", "ncls"):
            del sub["stabilization"]
    split = {s["name"]: s for s in split_problem(doc, "level_set_five")["subproblems"]}
    assert split["vof"]["stabilization"] == {"method": "SUPG", "tau": "element_size"}
    # a sub-problem that states one does not merge with the entry's: flow is ASGS with its own keys
    assert split["flow"]["stabilization"]["method"] == "ASGS" and split["flow"]["stabilization"]["tau"] == "metric"
    assert split["mcorr"]["stabilization"] == {"method": "SUPG", "tau": "element_size"}


# --- checked on load ---------------------------------------------------------------------------------------------

def _state(doc, sub, block):
    subproblem(doc, "level_set_five", sub)["stabilization"] = block
    return doc


def test_a_bad_stabilization_is_refused_on_load_and_says_where():
    doc = _state(load_ymf(FORMULATION), "vof", {"method": "SUPG", "shock_capturing": {"factor": 0.5, "beta": 3.0}})
    with pytest.raises(VocabularyError, match=r"subproblems\[vof\]\.stabilization\.shock_capturing\.beta: 3 is outside"):
        check(doc)


def test_the_strict_schema_refuses_a_key_it_does_not_know(tmp_path):
    text = FORMULATION.read_text(encoding="utf-8").replace("tau: metric\n            lagged: false\n            shock_capturing: {factor: 0.5, beta: 1.5",
                                                           "tau: metric\n            damping: 1.0\n            shock_capturing: {factor: 0.5, beta: 1.5", 1)
    (tmp_path / "dambreak.ymf").write_text((TWO_PHASE / "dambreak.ymf").read_text(encoding="utf-8"), encoding="utf-8")
    (tmp_path / "two_phase_level_set.ymf").write_text(
        (TWO_PHASE / "two_phase_level_set.ymf").read_text(encoding="utf-8"), encoding="utf-8")
    (tmp_path / "bad.ymf").write_text(text, encoding="utf-8")
    with pytest.raises(Exception, match="damping"):
        load_ymf(tmp_path / "bad.ymf")


def test_a_method_word_alone_still_loads():
    doc = load_ymf(FORMULATION)
    for sub in doc["solution_paths"]["discretizations"][0]["subproblems"]:
        if "stabilization" in sub:
            sub["stabilization"] = "SUPG"
    check(doc)


# --- the archive ---------------------------------------------------------------------------------------------------

def test_the_stabilization_survives_the_archive_and_is_part_of_its_key(tmp_path):
    realization = {"h": 0.0322}
    spec, again, digest = roundtrip(tmp_path, FORMULATION, "level_set_five", realization)
    for name in ("flow", "vof", "ncls", "rdls"):
        assert subproblem(again, "level_set_five", name)["stabilization"] == subproblem(
            spec, "level_set_five", name)["stabilization"]
    shock = subproblem(again, "level_set_five", "vof")["stabilization"]["shock_capturing"]
    assert shock == {"factor": 0.5, "beta": 1.5, "reference": 1.0, "lagged": True}
    # it is an input: another value is another key, and the spec that states none has another key than this one
    from ymf import closure
    other = load_ymf(FORMULATION)
    _state(other, "vof", {"method": "SUPG", "shock_capturing": {"factor": 0.6}})
    assert closure.input_digest(other, "level_set_five", realization) != digest
    plain, _, _ = closure.load(DAMBREAK)
    assert closure.input_digest(plain, "level_set_five", realization) != digest


def test_a_spec_in_the_older_spelling_has_the_words_of_the_newer(tmp_path):
    old = FORMULATION.read_text(encoding="utf-8").replace("lagged:", "lag:").replace("beta:", "sc_beta:").replace(
        "reference:", "sc_uref:")
    # the replacement also changed the (unused) keys of sub-problems that do not read them; the old spellings are accepted
    for name in ("dambreak.ymf", "two_phase_level_set.ymf"):
        (tmp_path / name).write_text((TWO_PHASE / name).read_text(encoding="utf-8"), encoding="utf-8")
    (tmp_path / "old.ymf").write_text(old, encoding="utf-8")
    new = split_problem(load_ymf(FORMULATION), "level_set_five")
    older = split_problem(load_ymf(tmp_path / "old.ymf"), "level_set_five")
    assert [s.get("stabilization") for s in older["subproblems"]] == [s.get("stabilization") for s in new["subproblems"]]


def test_rdls_stated_at_beta_1_without_a_reference_is_complete_and_round_trips_an_archive(tmp_path):
    from ymf import closure
    from ymf.symbolic import split_problem
    doc = load_ymf(FORMULATION)
    rdls = {"method": "SUPG", "tau": "metric", "lagged": False,
            "shock_capturing": {"factor": 0.9, "beta": 1.0, "lagged": False}}
    _state(doc, "rdls", rdls)
    check(doc)
    path = tmp_path / "rdls_beta1.archive.ymf"
    closure.write(path, doc, {})
    again, _, _ = closure.load(path)
    got = {s["name"]: s for s in split_problem(again, "level_set_five")["subproblems"]}["rdls"]["stabilization"]
    assert got == resolve(rdls) and "reference" not in got["shock_capturing"] and missing(got) == []
    assert closure.input_digest(again, "level_set_five", {"h": 0.0322}) == closure.input_digest(
        doc, "level_set_five", {"h": 0.0322})
