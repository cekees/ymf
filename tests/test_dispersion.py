"""ymf.dispersion (C42 round 5): the dispersive shallow-water family, solved by hyperbolic relaxation, stated so that C5's
dispatcher maps it onto GN_SW2DCV: the closure, the coefficients it supplies, and the three auxiliary unknowns."""

import copy

import pytest

pytest.importorskip("strictyaml")
pytest.importorskip("sympy")

from _formulation import SHALLOW_WATER, roundtrip, subproblem  # noqa: E402
from ymf.dispersion import AUXILIARY_ROLES, DispersionError, check, resolve  # noqa: E402
from ymf.references import YmfReferenceError  # noqa: E402
from ymf.schema import load_ymf  # noqa: E402
from ymf.state import reached  # noqa: E402
from ymf.symbolic import split_problem  # noqa: E402
from ymf.vocabulary import VocabularyError, check as check_vocabulary, coefficient_by_role, family_of  # noqa: E402

CIRCULAR = SHALLOW_WATER / "circular_dam_dispersive.ymf"
GN = "relaxed_convex_limited"

BLOCK = {"closure": "hyperbolic_relaxation", "closes": ["Pd", "Qd"], "auxiliaries": [
    {"name": "hη", "role": "relaxed_depth", "units": "m2"},
    {"name": "hw", "role": "relaxed_divergence"},
    {"name": "hβ", "role": "relaxed_bed_flux"}]}


def _block(**changes):
    out = copy.deepcopy(BLOCK)
    out.update(changes)
    return out


# --- resolve ------------------------------------------------------------------------------------------------------

def test_the_defaults_are_the_decks_and_each_auxiliary_says_whose_boundary_data_it_takes():
    r = resolve(BLOCK)
    assert [(a["name"], a["initial"], a["initial_is_default"], a["boundary_like"]) for a in r["auxiliaries"]] == [
        ("hη", "h^2", True, "h"), ("hw", "0", True, "h"), ("hβ", "0", True, "q")]
    assert r["closes"] == ["Pd", "Qd"] and r["closure"] == "hyperbolic_relaxation"


def test_a_stated_initial_value_replaces_the_default_and_says_so():
    # a solitary wave's hw is not 0 (test/SWFlow/santos_step.py)
    aux = copy.deepcopy(BLOCK["auxiliaries"])
    aux[1]["initial"] = "-h^2 c h0 dh"
    hw = resolve(_block(auxiliaries=aux))["auxiliaries"][1]
    assert (hw["initial"], hw["initial_is_default"]) == ("-h^2 c h0 dh", False)


def test_the_three_roles_are_the_ones_of_the_decks_comment():
    # "heta -> h^2, hbeta -> q.grad(Z), hw -> h^2 div(u) + 3/2 hbeta; the BCs for heta and hw are those of h, hbeta's those of x_mom"
    assert {r: (a["consistent"], a["boundary_like"]) for r, a in AUXILIARY_ROLES.items()} == {
        "relaxed_depth": ("h^2", "h"), "relaxed_divergence": ("h^2 ∇·u + (3/2) hβ", "h"),
        "relaxed_bed_flux": ("q·∇Z", "q")}


@pytest.mark.parametrize("block, message", [
    ("hyperbolic_relaxation", "give closure"),
    ({"closes": ["Pd"]}, "give closure"),
    (_block(closure="projection"), "closure: projection is not one of hyperbolic_relaxation"),
    (_block(extra=1), "no key extra"),
    (_block(closes=[]), "closes: names the model's coefficients"),
    (_block(closes=["Pd", "Pd"]), "closes: a coefficient twice"),
    (_block(auxiliaries=[]), "one per role"),
    (_block(auxiliaries=BLOCK["auxiliaries"][:2]), "missing relaxed_bed_flux"),
    (_block(auxiliaries=BLOCK["auxiliaries"] + [{"name": "hq", "role": "relaxed_depth"}]),
     r"auxiliaries\[3\]\.role: relaxed_depth is already the role of hη"),
    (_block(auxiliaries=[BLOCK["auxiliaries"][0], {"name": "hη", "role": "relaxed_divergence"},
                         BLOCK["auxiliaries"][2]]), "name: hη is already an auxiliary"),
    (_block(auxiliaries=[{"name": "hη", "role": "temperature"}] + BLOCK["auxiliaries"][1:]),
     "role: temperature is not one of relaxed_depth"),
    (_block(auxiliaries=[{"role": "relaxed_depth"}] + BLOCK["auxiliaries"][1:]), "give name: and role:"),
    (_block(auxiliaries=[dict(BLOCK["auxiliaries"][0], damping=1)] + BLOCK["auxiliaries"][1:]), "no key damping"),
])
def test_a_dispersion_that_does_not_make_sense_is_refused(block, message):
    with pytest.raises(DispersionError, match=message):
        resolve(block)


# --- the spec -----------------------------------------------------------------------------------------------------

@pytest.fixture(scope="module")
def circular():
    return load_ymf(CIRCULAR)


def test_the_model_is_the_dispersive_family_with_the_two_coefficients_the_closure_supplies(circular):
    assert family_of(circular) == "dispersive_shallow_water"
    assert circular["model"]["physical_model"]["assumptions"] == ["depth_averaged", "dispersive"]
    assert coefficient_by_role(circular, "dispersive_pressure") == "Pd"
    assert coefficient_by_role(circular, "dispersive_bed_force") == "Qd"
    # the family is abstract, the case reaches discretized with cells per axis not yet needed: one square count
    assert reached(circular)[0] == "discretized"


def test_the_split_problem_gives_the_six_unknowns_of_gn_sw2dcv(circular):
    (gn,) = split_problem(circular, GN)["subproblems"]
    # h and the two components of q from the model, and the three auxiliaries the closure introduces
    assert gn["unknowns"] == {"h": ["h"], "q": ["q_0", "q_1"]}
    aux = [a["name"] for a in gn["dispersion"]["auxiliaries"]]
    assert aux == ["hη", "hw", "hβ"]
    assert sum(len(v) for v in gn["unknowns"].values()) + len(aux) == 6
    # the closure supplies the dispersive pressure and bed force at each point, with their gradients (a force): the
    # momentum reads the gradient of Pd, the mass nothing
    assert gn["reads"]["Pd"] == {"via": "dispersion", "at": "current_step", "components": ["Pd"]}
    assert gn["reads"]["Qd"]["components"] == ["Qd"]
    momentum = {e["component"]: e for e in gn["adr"]["equations"]}
    assert "grad_Pd_0" in momentum["q_0"]["reaction"]["r"] and "grad_Pd_1" in momentum["q_1"]["reaction"]["r"]
    assert "reaction" not in momentum["h"]
    assert gn["coefficient_roles"]["dispersive_pressure"] == "Pd"


def test_the_circular_dam_states_what_the_deck_set_by_hand(circular):
    # test/SWFlow/circular_damBreak.py, sw_model 1: hη = h², hw = hβ = 0, walls all round, 49 × 49 cells
    (gn,) = split_problem(circular, GN)["subproblems"]
    assert {a["role"]: a["initial"] for a in gn["dispersion"]["auxiliaries"]} == {
        "relaxed_depth": "h^2", "relaxed_divergence": "0", "relaxed_bed_flux": "0"}
    assert all(c["variable"] == "q" and c["type"] in ("dirichlet", "neumann") for c in gn["boundary_conditions"])
    assert [e["mesh"]["cells"] for e in circular["solution_paths"]["discretizations"]] == [49]


# --- checked on load ----------------------------------------------------------------------------------------------

def _gn(doc):
    return subproblem(doc, GN, "gn")


def test_a_model_that_is_not_dispersive_has_no_dispersion(tmp_path):
    doc = load_ymf(SHALLOW_WATER / "torino.ymf")
    entry = doc["solution_paths"]["discretizations"][0]
    entry["subproblems"][0]["dispersion"] = copy.deepcopy(BLOCK)
    with pytest.raises(VocabularyError, match="only a model of family dispersive_shallow_water has a dispersion; this one is shallow_water"):
        check_vocabulary(doc)


@pytest.mark.parametrize("change, message", [
    (lambda d: d.__setitem__("closes", ["Pd", "g"]), "g has role gravity; hyperbolic_relaxation supplies"),
    (lambda d: d.__setitem__("closes", ["Pd"]), "supplies dispersive_pressure and dispersive_bed_force, and these are dispersive_pressure"),
    (lambda d: d.__setitem__("closes", ["Pd", "nothing"]), "nothing has role none"),
    (lambda d: d["auxiliaries"][0].__setitem__("name", "h"), "the auxiliary h is already a name in the model"),
    (lambda d: d["auxiliaries"][1].__setitem__("name", "Pd"), "the auxiliary Pd is already a name in the model"),
])
def test_a_dispersion_that_disagrees_with_the_model_is_refused_on_load_and_says_where(circular, change, message):
    doc = copy.deepcopy(circular)
    block = _gn(doc)["dispersion"]
    block["auxiliaries"] = [dict(a) for a in block["auxiliaries"]]
    change(block)
    with pytest.raises(VocabularyError, match=r"subproblems\[gn\]\.dispersion.*: .*" + message):
        check_vocabulary(doc)


def test_the_strict_schema_refuses_a_key_it_does_not_know(tmp_path):
    for name in ("shallow_water.ymf", "dispersive.ymf", "dispersive_numerics.ymf"):
        (tmp_path / name).write_text((SHALLOW_WATER / name).read_text(encoding="utf-8"), encoding="utf-8")
    text = (SHALLOW_WATER / "circular_dam_dispersive.ymf").read_text(encoding="utf-8")
    (tmp_path / "bad.ymf").write_text(text, encoding="utf-8")
    bad = (tmp_path / "dispersive_numerics.ymf")
    bad.write_text(bad.read_text(encoding="utf-8").replace("closes: [Pd, Qd]", "closes: [Pd, Qd]\n            relaxation: 1.0"),
                   encoding="utf-8")
    with pytest.raises(Exception, match="relaxation"):
        load_ymf(tmp_path / "bad.ymf")


# --- the archive --------------------------------------------------------------------------------------------------

def test_the_dispersion_survives_the_archive_and_is_part_of_its_key(tmp_path):
    from ymf import closure
    realization = {"cells": 49}
    spec, again, digest = roundtrip(tmp_path, CIRCULAR, GN, realization)
    assert resolve(_gn(again)["dispersion"]) == resolve(_gn(spec)["dispersion"])
    changed = load_ymf(CIRCULAR)
    _gn(changed)["dispersion"]["auxiliaries"][1]["initial"] = "0.5"
    assert closure.input_digest(changed, GN, realization) != digest
