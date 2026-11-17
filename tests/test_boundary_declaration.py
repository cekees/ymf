"""Boundary kind (strong, weak, none) and porosity (C42 round 5): how a condition reaches the discrete problem, and the
porosity of the two-phase model, stated so that a generator picks its face kernel and its coefficient by the spec."""

import pytest

pytest.importorskip("strictyaml")
pytest.importorskip("sympy")

from _formulation import TWO_PHASE, roundtrip, subproblem  # noqa: E402
from ymf.boundaries import BoundaryError, check_condition, kind_of, kinds  # noqa: E402
from ymf.conditions import expand_condition  # noqa: E402
from ymf.schema import load_ymf  # noqa: E402
from ymf.symbolic import split_problem  # noqa: E402
from ymf.vocabulary import VocabularyError, check, coefficient_by_role  # noqa: E402

DAMBREAK = TWO_PHASE / "dambreak.ymf"
FORMULATION = TWO_PHASE / "dambreak_formulation.ymf"
POROUS = TWO_PHASE / "dambreak_porous.ymf"


# --- the kind of a condition --------------------------------------------------------------------------------------

@pytest.mark.parametrize("condition, kind", [
    ({"type": "dirichlet"}, "strong"),
    ({"type": "dirichlet", "imposed": "weak"}, "weak"),
    ({"type": "slip"}, "strong"), ({"type": "inflow"}, "strong"), ({"type": "characteristic"}, "strong"),
    ({"type": "neumann"}, "none"), ({"type": "robin"}, "none"), ({"type": "open"}, "none"),
    ({"type": "open", "imposed": "none"}, "none"),
    ({"type": "periodic"}, None),
])
def test_the_kind_is_stated_or_the_types_default(condition, kind):
    check_condition(condition, "c")
    assert kind_of(condition) == kind


@pytest.mark.parametrize("condition, message", [
    ({"type": "neumann", "imposed": "strong"}, "a neumann condition is imposed none"),
    ({"type": "robin", "imposed": "weak"}, "a robin condition is imposed none"),
    ({"type": "open", "imposed": "weak"}, "a open condition is imposed none"),
    ({"type": "dirichlet", "imposed": "none"}, "imposes a value, strong or weak"),
    ({"type": "periodic", "imposed": "weak"}, "periodic condition identifies faces"),
    ({"type": "dirichlet", "penalty": 10.0}, "penalty belong to a weak condition"),
    ({"type": "dirichlet", "imposed": "strong", "form": "symmetric"}, "form belong to a weak condition"),
    ({"type": "dirichlet", "imposed": "weak", "penalty": -1.0}, "penalty -1.0 is not a positive number"),
    ({"type": "dirichlet", "imposed": "weak", "form": "complete"}, "form complete is not one of symmetric"),
    ({"type": "dirichlet", "imposed": "maybe"}, "imposed maybe is not one of strong, weak, none"),
])
def test_a_kind_that_does_not_fit_its_type_is_refused(condition, message):
    with pytest.raises(BoundaryError, match=message):
        check_condition(condition, "c")


def test_an_alias_gives_its_kind_to_the_dirichlet_part_only():
    slip = expand_condition({"variable": "v", "type": "slip", "region": "w", "imposed": "weak", "penalty": 50.0,
                             "form": "incomplete"}, rank=1)
    normal, tangential = slip
    assert normal["type"] == "dirichlet" and (normal["imposed"], normal["penalty"], normal["form"]) == (
        "weak", 50.0, "incomplete")
    assert tangential["type"] == "neumann" and not {"imposed", "penalty", "form"} & set(tangential)
    assert kind_of(tangential) == "none"


# --- the dam break ------------------------------------------------------------------------------------------------

def test_the_dam_break_states_how_each_condition_is_imposed():
    by = {(k["region"], k["variable"], k["type"]): k for k in kinds(load_ymf(FORMULATION))}
    assert by[("walls", "v", "dirichlet")]["kind"] == "weak" and by[("walls", "v", "dirichlet")]["form"] == "symmetric"
    assert by[("walls", "v", "neumann")]["kind"] == "none"
    assert by[("top", "p", "dirichlet")]["kind"] == "weak"
    assert by[("top", "v", "neumann")]["kind"] == "none"
    assert by[("top", "χ", "dirichlet")]["kind"] == "strong" and by[("top", "χ", "dirichlet")]["only"] == "inflow"


def test_the_split_problem_gives_each_sub_problem_its_kinds_so_a_generator_picks_its_faces():
    split = {s["name"]: s for s in split_problem(load_ymf(FORMULATION), "level_set_five")["subproblems"]}
    # VOF: the air entering through the top is strong, an upwind flux with the Dirichlet inflow (the engine's Strong)
    (inflow,) = split["vof"]["boundary_conditions"]
    assert (inflow["kind"], inflow["only"], inflow["carried_by"]) == ("strong", "inflow", "χh")
    # RANS2P: weak Dirichlet with the symmetric adjoint (WeakDirichlet), the open part natural (Flux)
    assert {(c["variable"], c["type"], c["kind"]) for c in split["flow"]["boundary_conditions"]} == {
        ("v", "dirichlet", "weak"), ("v", "neumann", "none"), ("p", "dirichlet", "weak")}
    assert all(c["form"] == "symmetric" for c in split["flow"]["boundary_conditions"] if c["kind"] == "weak")


def test_a_spec_that_states_no_kind_gets_the_types_default_in_the_split_problem():
    split = {s["name"]: s for s in split_problem(load_ymf(DAMBREAK), "level_set_five")["subproblems"]}
    assert {c["kind"] for c in split["flow"]["boundary_conditions"]} == {"strong", "none"}
    assert not any("form" in c or "penalty" in c for c in split["flow"]["boundary_conditions"])


def test_a_bad_kind_in_a_condition_set_is_refused_on_load_and_says_where():
    doc = load_ymf(FORMULATION)
    doc["model"]["strong_form"]["condition_sets"]["atmosphere"][1] = {"variable": "v", "type": "open", "imposed": "weak"}
    with pytest.raises(VocabularyError, match=r"condition_sets\.atmosphere\[1\] \(open on v\).*imposed none"):
        check(doc)


def test_the_kinds_survive_the_archive_and_are_part_of_its_key(tmp_path):
    from ymf import closure
    realization = {"h": 0.0322}
    spec, again, digest = roundtrip(tmp_path, FORMULATION, "level_set_five", realization)
    assert kinds(again) == kinds(spec)
    changed = load_ymf(FORMULATION)
    changed["model"]["strong_form"]["condition_sets"]["atmosphere"][2]["imposed"] = "weak"
    assert closure.input_digest(changed, "level_set_five", realization) != digest


# --- porosity -----------------------------------------------------------------------------------------------------

@pytest.fixture(scope="module")
def porous():
    return load_ymf(POROUS)


def test_the_two_phase_model_states_its_porosity_by_role(porous):
    assert coefficient_by_role(porous, "porosity") == "θ"
    assert coefficient_by_role(load_ymf(DAMBREAK), "porosity") is None


def test_the_sub_problems_say_which_use_the_porosity(porous):
    split = {s["name"]: s for s in split_problem(porous, "level_set_five")["subproblems"]}
    assert {n for n, s in split.items() if "porosity" in s.get("coefficient_roles", {})} == {"flow", "vof", "mcorr"}
    assert split["flow"]["coefficient_roles"]["density"] == "ρ"
    # NCLS has none, and the plain dam break none anywhere
    plain = {s["name"]: s for s in split_problem(load_ymf(DAMBREAK), "level_set_five")["subproblems"]}
    assert not any("porosity" in s.get("coefficient_roles", {}) for s in plain.values())


def test_vof_carries_the_porosity_in_its_mass_and_its_flux_as_the_generator_spikes_hand_made_form_did(porous):
    # generator/make_inputs.py wrote "∂(ε*χ)/∂t + ∇·(ε*χ*v) = 0" in a Space by hand (porosity ε a per-point field);
    # the spec now states it, with θ inlined as the expression of the coordinates it is
    split = {s["name"]: s for s in split_problem(porous, "level_set_five")["subproblems"]}
    eq = split["vof"]["adr"]["equations"][0]
    assert eq["mass"]["m"].startswith("χh*(") and "numpy.select" in eq["mass"]["m"]
    numpy = pytest.importorskip("numpy")
    env = {"numpy": numpy, "h": 0.0322, "x": 2.3, "y": 0.2, "χh": 1.0}
    run = lambda code: eval(code, {"__builtins__": {}}, env)  # noqa: S307
    assert run(eq["mass"]["dm"]["χh"]) == pytest.approx(run(eq["mass"]["m"]))      # m = θ χh, so dm/dχh = θ
    assert all(f.startswith("χh*") and "v_%d" % i in f for i, f in enumerate(eq["advection"]["f"]))
    assert eq["mass"]["depends_on"] == {"χh": "linear"}
    # NCLS is the dam break's, unchanged
    plain = {s["name"]: s for s in split_problem(load_ymf(DAMBREAK), "level_set_five")["subproblems"]}
    assert split["ncls"]["adr"] == plain["ncls"]["adr"]


def test_the_porosity_is_one_where_it_is_one_and_below_one_in_the_block(porous):
    numpy = pytest.importorskip("numpy")
    split = {s["name"]: s for s in split_problem(porous, "level_set_five")["subproblems"]}
    m = split["vof"]["adr"]["equations"][0]["mass"]["dm"]["χh"]
    env = {"numpy": numpy, "h": 0.0322}
    at = lambda x, y: eval(m, {"__builtins__": {}}, dict(env, x=x, y=y))  # noqa: S307
    assert at(1.0, 1.0) == pytest.approx(1.0)          # open flow
    assert at(3.0, 0.1) == pytest.approx(0.6)          # inside the block
    assert 0.6 < at(2.2, 0.1) < 1.0                    # the edge, smoothed over 1.5 h


def test_the_porosity_survives_the_archive(tmp_path):
    spec, again, _ = roundtrip(tmp_path, POROUS, "level_set_five", {"h": 0.0322})
    assert coefficient_by_role(again, "porosity") == "θ"
    assert subproblem(again, "level_set_five", "vof")["equation"] == "∂(θ χh)/∂t + ∇·(θ v χh) = 0"
