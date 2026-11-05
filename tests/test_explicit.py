"""A closed spec writes out what the input leaves to rules (C42 round 2).

ymf.explicit: the coefficients of the state named with the unknowns they
follow (depends_on), and every sub-problem of a split discretization with
its unknown (for), solver, regularization and lags in full (carried_by
for a coefficient read through its carrier). Like the boundary-condition
expansion (#9), it is idempotent, leaves a spec without shorthand alone,
and keys a spec and the archive it closes into alike.
"""

import copy
import glob
from pathlib import Path

import pytest

pytest.importorskip("strictyaml")

from ymf import closure  # noqa: E402
from ymf.compose import load_composed  # noqa: E402
from ymf.explicit import explicit, state_coefficients, without_derived  # noqa: E402
from ymf.names import FUNCTION_NAMES, names_in  # noqa: E402
from ymf.references import YmfReferenceError, resolve  # noqa: E402

ROOT = Path(__file__).parent.parent
FIXTURES = Path(__file__).parent / "fixtures"
DAMBREAK = FIXTURES / "two_phase" / "dambreak.ymf"
TORINO = FIXTURES / "shallow_water" / "torino.ymf"


def _documents():
    for path in sorted(glob.glob(str(ROOT / "examples" / "**" / "*.ymf"), recursive=True)
                       + glob.glob(str(FIXTURES / "**" / "*.ymf"), recursive=True)):
        try:
            yield Path(path), load_composed(path, check_state=False)
        except Exception:          # an input made to fail, or a fragment
            continue


def _subproblem(doc, name, discretization="level_set_five"):
    entry = next(d for d in doc["solution_paths"]["discretizations"]
                 if d["name"] == discretization)
    return next(s for s in entry["subproblems"] if s["name"] == name)


# -- the name reader -----------------------------------------------------------

def test_names_in_reads_a_formula_as_the_parser_does():
    declared = {"ρw", "ρa", "χ", "φ", "h", "q", "n", "g", "x", "y", "t"}
    assert names_in("ρw (1 − χ) + ρa χ", declared) == {"ρw", "ρa", "χ"}
    assert names_in("heaviside(φ)", declared) == {"φ"}
    # a function call is not a product of letters (h, n in norm and heaviside)
    assert names_in("g n² norm(q) q / h^(7/3)", declared) == {"g", "n", "q", "h"}
    # an undeclared run of letters is the product of its letters
    assert names_in("gh", declared) == {"g", "h"}


def test_the_name_readers_functions_are_the_notations():
    sympy = pytest.importorskip("sympy")  # noqa: F841
    from ymf.symbolic.language import Space
    from ymf.symbolic.notation import Parser
    assert set(Parser(Space(3, {})).functions) == set(FUNCTION_NAMES)


def test_state_coefficients_agree_with_the_symbolic_reader_where_it_resolves():
    pytest.importorskip("sympy")
    from ymf.symbolic.problem import state_coefficients as symbolic
    checked = 0
    for path, doc in _documents():
        strong = doc["model"]["strong_form"]
        unknowns = [u if isinstance(u, str) else u["name"] for u in strong["unknowns"]]
        try:
            resolved = symbolic(strong)
        except Exception:
            continue
        text = state_coefficients(doc)
        # the symbolic reader leaves out what does not resolve (an abstract
        # model's ρ without ρw); where it resolves, the two agree
        for name, components in resolved.items():
            assert text[name] == sorted({c.rsplit("_", 1)[0] if c not in unknowns else c
                                         for c in components}, key=unknowns.index), path
            checked += 1
    assert checked >= 9


def test_an_abstract_model_names_its_state_coefficients_without_values():
    doc = load_composed(FIXTURES / "two_phase" / "two_phase_level_set.ymf")
    assert state_coefficients(doc) == {"χ": ["φ"], "ρ": ["φ"], "μ": ["φ"]}


# -- explicit() ------------------------------------------------------------------

def test_explicit_is_idempotent_and_leaves_plain_specs_alone():
    plain, written = 0, 0
    for path, doc in _documents():
        once = explicit(doc)
        assert explicit(once) == once, path
        if not state_coefficients(doc) and not any(
                d.get("subproblems") for d in
                (doc.get("solution_paths") or {}).get("discretizations") or []):
            from ymf.conditions import expand
            assert once == expand(doc), path
            plain += 1
        else:
            written += 1
    assert plain >= 8 and written >= 4


def test_the_closed_dam_break_names_its_coefficients_of_the_state():
    closed = explicit(load_composed(DAMBREAK))
    coefficients = closed["model"]["strong_form"]["coefficients"]
    assert coefficients["χ"] == {"formula": "heaviside(φ)", "units": "1", "depends_on": ["φ"]}
    assert coefficients["ρ"]["depends_on"] == ["φ"]
    assert coefficients["μ"]["depends_on"] == ["φ"]
    assert "depends_on" not in coefficients["ρw"]


def test_the_closed_dam_break_writes_its_split_out():
    spec = load_composed(DAMBREAK)
    closed = explicit(spec)
    flow = _subproblem(closed, "flow")
    assert flow["lags"] == {"φ": {"at": "previous_step"},
                            "χ": {"via": "vof", "at": "previous_step", "carried_by": "χh"}}
    # the entry's solver, refined by the sub-problem's own
    assert flow["solver"] == {"type": "nonlinear", "nonlinear_solver": "newton",
                              "linear_solver": "petsc", "max_iterations": 100}
    assert flow["regularize"] == {"heaviside": {"by": "smoothed_heaviside", "width": "ε"}}
    vof = _subproblem(closed, "vof")
    assert vof["for"] == "χh"
    assert vof["solver"]["max_iterations"] == 50
    # redistancing's own sign regularization joins the entry's heaviside
    rdls = _subproblem(closed, "rdls")
    assert set(rdls["regularize"]) == {"heaviside", "sign"}
    assert rdls["regularize"]["sign"]["width"] == "εr"
    assert rdls["lags"]["φs"] == {"of": "φ", "via": "ncls", "at": "current_step"}
    assert _subproblem(closed, "mcorr")["for"] == "φc"
    # the input is not changed
    assert _subproblem(spec, "flow")["lags"]["φ"] == "previous_step"


def test_the_closed_shallow_water_case_names_friction_and_drag():
    closed = explicit(load_composed(TORINO))
    coefficients = closed["model"]["strong_form"]["coefficients"]
    for name in ("τb", "Fd", "Fr", "β"):
        assert coefficients[name]["depends_on"] == ["h", "q"], name
    assert "depends_on" not in coefficients["θ"]
    swe = _subproblem(closed, "swe", "convex_limited")
    assert swe["solver"] == {"type": "time_marching", "linear_solver": "petsc"}


@pytest.mark.parametrize("path", [DAMBREAK, TORINO, FIXTURES / "shallow_water" / "lake_at_rest.ymf",
                                  ROOT / "examples" / "nonlinear_diffusion.ymf"])
def test_a_closed_spec_reloads_with_its_specs_keys(tmp_path, path):
    spec, outputs, _ = closure.load(path, check_state=False)
    closure.write(tmp_path / "closed.ymf", spec, outputs)
    again, _, _ = closure.load(tmp_path / "closed.ymf", check_state=False)
    assert [k for k, _, _ in closure.planned_outputs(again)] == \
        [k for k, _, _ in closure.planned_outputs(spec)]
    assert explicit(again)["solution_paths"] == explicit(spec)["solution_paths"]
    assert explicit(again)["model"]["strong_form"]["coefficients"] == \
        explicit(spec)["model"]["strong_form"]["coefficients"]


def test_depends_on_does_not_count_towards_a_key():
    # nonlinear_diffusion's keys before C42 round 2 (cekees/ymf main b77262a)
    spec, _, _ = closure.load(ROOT / "examples" / "nonlinear_diffusion.ymf")
    keys = [k for k, _, _ in closure.planned_outputs(spec)]
    assert keys[:3] == ["P1/cells=8/55bb9712a3fc", "P1/cells=16/2b0926f8753a",
                        "P1/cells=32/f7ada584d9ca"]
    assert without_derived(explicit(spec))["model"] == spec["model"]


def test_split_problem_reads_the_closed_spec_as_the_spec():
    pytest.importorskip("sympy")
    from ymf.symbolic import split_problem
    spec = load_composed(DAMBREAK)
    assert split_problem(explicit(spec), "level_set_five") == \
        split_problem(spec, "level_set_five")


# -- what load refuses ------------------------------------------------------------

def test_a_false_depends_on_is_refused():
    doc = copy.deepcopy(load_composed(DAMBREAK))
    doc["model"]["strong_form"]["coefficients"]["ρ"]["depends_on"] = ["v"]
    with pytest.raises(YmfReferenceError, match="ρ.depends_on says"):
        resolve(doc)
    doc["model"]["strong_form"]["coefficients"]["ρ"]["depends_on"] = ["φ"]
    resolve(doc)


def test_a_false_carried_by_is_refused():
    doc = copy.deepcopy(load_composed(DAMBREAK))
    _subproblem(doc, "flow")["lags"]["χ"]["carried_by"] = "φc"
    with pytest.raises(YmfReferenceError, match="χh carries it"):
        resolve(doc)
    _subproblem(doc, "flow")["lags"]["φ"] = {"at": "previous_step", "carried_by": "χh"}
    _subproblem(doc, "flow")["lags"]["χ"]["carried_by"] = "χh"
    with pytest.raises(YmfReferenceError, match="no sub-problem introduces"):
        resolve(doc)
