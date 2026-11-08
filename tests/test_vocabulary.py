"""ymf's vocabulary (C42 round 3): the words a spec states and a solver's
dispatcher maps -- model families, sub-problem roles, the discretization's
words -- checked on load, derived ones written by closure.
"""

import copy
from pathlib import Path

import pytest

pytest.importorskip("strictyaml")

from ymf import vocabulary as v  # noqa: E402
from ymf.closure import input_digest  # noqa: E402
from ymf.explicit import explicit, without_derived  # noqa: E402
from ymf.references import YmfReferenceError, resolve  # noqa: E402
from ymf.schema import load_ymf  # noqa: E402

FIXTURES = Path(__file__).parent / "fixtures"
TWO_PHASE = FIXTURES / "two_phase"
SHALLOW = FIXTURES / "shallow_water"


@pytest.fixture
def dambreak():
    return load_ymf(TWO_PHASE / "dambreak.ymf")


@pytest.fixture
def clsvof():
    return load_ymf(TWO_PHASE / "dambreak_clsvof.ymf")


def _entry(doc, name):
    return next(e for e in doc["solution_paths"]["discretizations"] if e["name"] == name)


def _sub(doc, entry, name):
    return next(s for s in _entry(doc, entry)["subproblems"] if s["name"] == name)


# --- families ----------------------------------------------------------------

def test_the_families_are_identified_by_processes_and_assumptions(dambreak):
    assert v.family_of(dambreak) == "two_phase_level_set"
    assert v.family_of(load_ymf(SHALLOW / "lake_at_rest.ymf")) == "shallow_water"
    assert v.family_of(load_ymf(Path(__file__).parents[1] / "examples" / "poisson.ymf")) is None


def test_the_dispersive_family_and_its_older_key():
    doc = {"model": {"physical_model": {"assumptions": ["depth_averaged", "dispersive"]}}}
    assert v.family_of(doc) == "dispersive_shallow_water"
    # the dispatch table keyed it on how GN_SW2DCV discretizes it
    old = {"model": {"physical_model": {"assumptions": ["depth_averaged", "hyperbolic_relaxation"]}}}
    assert v.family_of(old) == "dispersive_shallow_water"


def test_a_stated_family_is_checked(dambreak):
    doc = copy.deepcopy(dambreak)
    doc["model"]["physical_model"]["family"] = "two_phase_level_set"
    resolve(doc)
    doc["model"]["physical_model"]["family"] = "shallow_water"
    with pytest.raises(YmfReferenceError, match="family is shallow_water, but its processes"):
        resolve(doc)
    doc["model"]["physical_model"]["family"] = "two_phase"
    with pytest.raises(YmfReferenceError, match="family two_phase: ymf knows"):
        resolve(doc)


# --- roles ---------------------------------------------------------------------

def test_the_five_model_split_has_the_dispatch_tables_roles(dambreak):
    assert v.roles_of(dambreak, _entry(dambreak, "level_set_five")) == {
        "flow": "velocity_pressure", "vof": "volume_fraction", "ncls": "level_set",
        "rdls": "redistancing", "mcorr": "mass_correction"}


def test_clsvof_solves_and_redistances_in_one_role(clsvof):
    assert v.roles_of(clsvof, _entry(clsvof, "level_set_clsvof")) == {
        "flow": "velocity_pressure", "clsvof": "level_set_redistanced"}


def test_the_projection_schemes_roles(dambreak):
    """damBreak.py's clsvof_projection: CLSVOF, then the velocity step, the
    pressure increment, the pressure update and a pressure spin-up."""
    subs = [{"name": "clsvof", "solves": ["interface"], "enforces": "norm(∇φ) = 1"},
            {"name": "flow", "solves": ["momentum"]},
            {"name": "pressureInc", "enforces": "∇·v = 0"},
            {"name": "pressure", "solves": ["continuity"]},
            {"name": "pressureInit", "spin_up": True}]
    assert [v.role_of(s, dambreak) for s in subs] == [
        "level_set_redistanced", "velocity", "pressure_increment", "pressure", "pressure_spin_up"]


def test_the_shallow_water_role():
    doc = load_ymf(SHALLOW / "lake_at_rest.ymf")
    assert v.roles_of(doc, doc["solution_paths"]["discretizations"][0]) == {"swe": "shallow_water"}


def test_a_stated_role_is_checked_against_the_structure(dambreak):
    doc = copy.deepcopy(dambreak)
    _sub(doc, "level_set_five", "vof")["role"] = "volume_fraction"
    _sub(doc, "level_set_five", "mcorr")["role"] = "mass_correction"
    resolve(doc)
    _sub(doc, "level_set_five", "ncls")["role"] = "redistancing"
    with pytest.raises(YmfReferenceError, match="ncls.*role says redistancing, but what it .* "
                                                "makes it level_set"):
        resolve(doc)
    _sub(doc, "level_set_five", "ncls")["role"] = "shallow_water"
    with pytest.raises(YmfReferenceError, match="not a role in family two_phase_level_set"):
        resolve(doc)
    # the old word for the spin-up is read as ymf's
    assert v.canonical("role", "pressure_initial") == "pressure_spin_up"


def test_closure_writes_the_family_and_the_roles_out_without_rekeying(dambreak):
    closed = explicit(dambreak)
    assert closed["model"]["physical_model"]["family"] == "two_phase_level_set"
    assert [s["role"] for s in _entry(closed, "level_set_five")["subproblems"]] == [
        "velocity_pressure", "volume_fraction", "level_set", "redistancing", "mass_correction"]
    assert explicit(closed) == closed                      # idempotent
    # derived: neither counts towards an output's key, stated or not
    assert without_derived(closed)["model"]["physical_model"].get("family") is None
    stated = copy.deepcopy(dambreak)
    stated["model"]["physical_model"]["family"] = "two_phase_level_set"
    _sub(stated, "level_set_five", "ncls")["role"] = "level_set"
    realization = {"h": 0.0322}
    assert (input_digest(stated, "level_set_five", realization)
            == input_digest(dambreak, "level_set_five", realization))


def test_a_model_of_no_family_closes_as_before():
    doc = load_ymf(Path(__file__).parents[1] / "examples" / "poisson.ymf")
    closed = explicit(doc)
    assert "family" not in (closed["model"].get("physical_model") or {})


# --- the discretization's words ------------------------------------------------

#: every label proteus's DispatchTable (cekees/proteus main 74d453d3) matches
#: on, with ymf's word for it: the table maps these, and C5 step 3 switches it
#: to the right-hand side
DISPATCH_TABLE_LABELS = [
    ("family", "two_phase_level_set", "two_phase_level_set"),
    ("family", "shallow_water", "shallow_water"),
    ("family", "dispersive_shallow_water", "dispersive_shallow_water"),
    ("role", "velocity_pressure", "velocity_pressure"),
    ("role", "velocity", "velocity"),
    ("role", "pressure_increment", "pressure_increment"),
    ("role", "pressure", "pressure"),
    ("role", "pressure_initial", "pressure_spin_up"),
    ("role", "volume_fraction", "volume_fraction"),
    ("role", "level_set", "level_set"),
    ("role", "level_set_redistanced", "level_set_redistanced"),
    ("role", "redistancing", "redistancing"),
    ("role", "mass_correction", "mass_correction"),
    ("role", "shallow_water", "shallow_water"),
    ("time_integration", "backward_euler", "backward_euler"),
    ("time_integration", "vbdf", "bdf2"),
    ("time_integration", "none", "none"),
    ("time_integration", "SSP_RK3", "ssp_rk3"),
    ("stabilization", "ASGS", "ASGS"),
    ("stabilization", "SUPG", "SUPG"),
    ("stabilization", "graph_viscosity_smoothness", "graph_viscosity_smoothness"),
    ("stabilization", "graph_viscosity_entropy_residual", "graph_viscosity_entropy_residual"),
    ("stabilization", "entropy_viscosity", "entropy_viscosity"),
    ("stabilization", "shock_capturing", "residual_shock_capturing"),
    ("scheme", "graph_viscosity_convex_limiting", "graph_viscosity_convex_limiting"),
    ("scheme_option", "entropy_viscosity_cE", "entropy_viscosity_coefficient"),
    ("scheme_option", "lumped_mass", "lumped_mass"),
    ("scheme_option", "bounds_relaxation", "bounds_relaxation"),
    ("step_control", "min_cfl", "min_cfl"),
    ("regularization", "smoothed_heaviside", "smoothed_heaviside"),
    ("binding_meaning", "element_size", "element_size"),
]


@pytest.mark.parametrize("kind, label, word", DISPATCH_TABLE_LABELS)
def test_every_dispatch_table_label_is_ymf_vocabulary(kind, label, word):
    assert v.canonical(kind, label) == word
    assert v.canonical(kind, word) == word


def test_the_inertia_forms_are_named():
    assert set(v.INERTIA_FORMS) == {"conservative_inertia", "nonconservative_inertia"}


@pytest.mark.parametrize("path, value, message", [
    (("subproblems", "flow", "time_integration"), "crank_nicolson",
     "flow].time_integration: time integration crank_nicolson: ymf knows backward_euler"),
    (("subproblems", "flow", "stabilization"), {"method": "GLS"},
     "flow].stabilization.method: stabilization GLS"),
    (("subproblems", "rdls", "writes_back"), {"φ": "overwrite"}, "writes_back.φ: writes back overwrite"),
    (("subproblems", "rdls", "regularize"), {"sign": {"by": "tanh", "width": "εr"}},
     "regularization tanh"),
    (("step_control",), {"type": "max_cfl"}, "step control max_cfl"),
    (("binds",), {"h": {"meaning": "cell_size"}, "ε": {"formula": "1.5 h"}}, "binding meaning cell_size"),
])
def test_a_word_ymf_does_not_have_is_refused(dambreak, path, value, message):
    doc = copy.deepcopy(dambreak)
    entry = _entry(doc, "level_set_five")
    if path[0] == "subproblems":
        _sub(doc, "level_set_five", path[1])[path[2]] = value
    else:
        entry[path[0]] = value
    with pytest.raises(YmfReferenceError, match=message):
        resolve(doc)


def test_the_old_words_still_load(dambreak):
    doc = copy.deepcopy(dambreak)
    _sub(doc, "level_set_five", "flow")["time_integration"] = "vbdf"
    _sub(doc, "level_set_five", "vof")["stabilization"] = {"method": "shock_capturing"}
    resolve(doc)


def test_a_scheme_and_its_options():
    doc = load_ymf(SHALLOW / "lake_at_rest.ymf")
    block = doc["solution_paths"]["discretizations"][0]["discretization"]
    assert block["scheme"] == "graph_viscosity_convex_limiting"
    resolve(doc)
    block["limiter"] = "minmod"
    with pytest.raises(YmfReferenceError, match="has no option limiter; its options are "
                                                "entropy_viscosity_coefficient"):
        resolve(doc)
    del block["limiter"]
    block["entropy_viscosity_coefficient"] = block.pop("entropy_viscosity_cE")
    resolve(doc)
    block["scheme"] = "muscl"
    with pytest.raises(YmfReferenceError, match="scheme muscl: ymf knows"):
        resolve(doc)
