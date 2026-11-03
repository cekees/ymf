"""A discretization split into sub-problems (C42).

tests/fixtures/two_phase holds the two-phase level-set model and the dam
break split the way proteus solves it: flow, VOF, level set, redistancing
and mass correction. ymf checks the split's structure and names on load; a
solver's dispatcher is what maps each sub-problem onto a module.
"""

import copy
from pathlib import Path

import pytest

pytest.importorskip("strictyaml")

from ymf.discretizations import subproblem_solver  # noqa: E402
from ymf.references import YmfReferenceError, resolve  # noqa: E402
from ymf.schema import load_ymf  # noqa: E402
from ymf.state import reached  # noqa: E402

TWO_PHASE = Path(__file__).parent / "fixtures" / "two_phase"


@pytest.fixture
def dambreak():
    return load_ymf(TWO_PHASE / "dambreak.ymf")


def _split(doc):
    return doc["solution_paths"]["discretizations"][0]


def _sub(doc, name):
    return next(s for s in _split(doc)["subproblems"] if s["name"] == name)


def test_the_model_alone_is_abstract():
    model = load_ymf(TWO_PHASE / "two_phase_level_set.ymf")
    assert model["ymf_state"] == "abstract"


def test_the_split_dam_break_is_discretized(dambreak):
    state, _ = reached(dambreak)
    assert state == "discretized"
    assert [s["name"] for s in _split(dambreak)["subproblems"]] == [
        "flow", "vof", "ncls", "rdls", "mcorr"]


def test_a_sub_problem_refines_the_entry_solver(dambreak):
    entry = _split(dambreak)
    assert subproblem_solver(entry, _sub(dambreak, "flow"))["max_iterations"] == 100
    assert subproblem_solver(entry, _sub(dambreak, "ncls"))["max_iterations"] == 50
    assert subproblem_solver(entry, _sub(dambreak, "ncls"))["nonlinear_solver"] == "newton"


def test_petsc_options_merge_key_by_key():
    entry = {"solver": {"linear_solver": "petsc",
                        "petsc_options": {"ksp_type": "gmres", "pc_type": "asm"}}}
    sub = {"solver": {"petsc_options": {"pc_type": "lu"}}}
    assert subproblem_solver(entry, sub)["petsc_options"] == {"ksp_type": "gmres",
                                                              "pc_type": "lu"}


def test_the_split_is_part_of_the_output_key(dambreak):
    from ymf.closure import planned_outputs
    before = [k for k, _, _ in planned_outputs(dambreak)]
    _sub(dambreak, "flow")["solver"]["max_iterations"] = 99
    after = [k for k, _, _ in planned_outputs(dambreak)]
    assert before != after


def test_the_whole_model_is_not_one_adr_system(dambreak):
    """v·∇φ is a Hamiltonian in another unknown: classified only per sub-problem."""
    pytest.importorskip("sympy")
    from ymf.symbolic import SymbolicError, adr_problem
    with pytest.raises(SymbolicError):
        adr_problem(dambreak)


def _broken(doc, change):
    doc = copy.deepcopy(doc)
    change(doc)
    return doc


def _set(path, value):
    def change(doc):
        node = doc
        for key in path[:-1]:
            node = node[key] if not isinstance(key, tuple) else next(
                s for s in node if s["name"] == key[0])
        node[path[-1]] = value
    return change


SPLIT = ("solution_paths", "discretizations", ("level_set_five",))


@pytest.mark.parametrize("change, message", [
    (_set(SPLIT + ("subproblems", ("ncls",), "solves"), ["interfaces"]),
     "names 'interfaces', which is not a model equation"),
    (_set(SPLIT + ("subproblems", ("ncls",), "solves"), ["momentum"]),
     "momentum is solved by both flow and ncls"),
    (_set(SPLIT + ("subproblems", ("ncls",), "solves"), []),
     "neither solves a model equation"),
    (_set(SPLIT + ("subproblems", ("rdls",), "name"), "flow"),
     "flow named more than once"),
    (_set(SPLIT + ("subproblems", ("vof",), "introduces"), {"name": "ρ"}),
     "ρ is already a name \\(in the model\\)"),
    (_set(SPLIT + ("subproblems", ("mcorr",), "introduces"), {"name": "φc",
                                                              "correction_of": "χ"}),
     "names 'χ', which is not a model unknown"),
    (_set(SPLIT + ("subproblems", ("ncls",), "lags"), {"w": "current_step"}),
     "names 'w', which is not a model unknown, coefficient or introduced unknown"),
    (_set(SPLIT + ("subproblems", ("ncls",), "lags"), {"v": {"via": "rans", "at": "current_step"}}),
     "names 'rans', which is not a sub-problem"),
    (_set(SPLIT + ("subproblems", ("rdls",), "writes_back"), {"ρ": "replace"}),
     "names 'ρ', which is not a model or introduced unknown"),
    (_set(SPLIT + ("binds",), {"g": {"value": 1.0}}),
     "g is already a name in the model"),
    (_set(SPLIT + ("regularize",), {"sqrt": {"by": "smooth"}}),
     "only the sharp functions heaviside, sign are regularized"),
])
def test_a_split_that_does_not_add_up_is_refused(dambreak, change, message):
    with pytest.raises(YmfReferenceError, match=message):
        resolve(_broken(dambreak, change))


def test_every_model_equation_is_solved_by_some_sub_problem(dambreak):
    doc = _broken(dambreak, lambda d: _split(d)["subproblems"].remove(_sub(d, "ncls")))
    with pytest.raises(YmfReferenceError, match="no sub-problem solves interface"):
        resolve(doc)


def test_splitting_needs_named_equations(dambreak):
    def unname(doc):
        strong = doc["model"]["strong_form"]
        strong["equations"] = [e["equation"] for e in strong["equations"]]
    doc = _broken(dambreak, unname)
    with pytest.raises(YmfReferenceError, match="name them"):
        resolve(doc)


def test_the_file_level_errors_reach_the_loader(tmp_path):
    text = (TWO_PHASE / "dambreak.ymf").read_text(encoding="utf-8")
    (tmp_path / "two_phase_level_set.ymf").write_text(
        (TWO_PHASE / "two_phase_level_set.ymf").read_text(encoding="utf-8"), encoding="utf-8")
    (tmp_path / "dambreak.ymf").write_text(
        text.replace("solves: [interface]", "solves: [interfaces]"), encoding="utf-8")
    with pytest.raises(YmfReferenceError, match="not a model equation"):
        load_ymf(tmp_path / "dambreak.ymf")
