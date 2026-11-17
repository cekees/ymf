"""Cells per axis (C42 round 5): ``mesh: {cells: {x: 80, y: 40}}``, so a structured mesh of a box that is not square is
stated, and a list stays a study."""

import copy

import pytest

pytest.importorskip("strictyaml")
pytest.importorskip("sympy")

from _formulation import FIXTURES, SHALLOW_WATER, roundtrip  # noqa: E402
from ymf import closure  # noqa: E402
from ymf.discretizations import CellsError, cells_per_axis, cells_text, domain_dim  # noqa: E402
from ymf.schema import load_ymf  # noqa: E402
from ymf.state import reached  # noqa: E402
from ymf.vocabulary import VocabularyError, check  # noqa: E402

DAM3BUMPS = SHALLOW_WATER / "dam3bumps.ymf"
BOX = "Ω = [0, 75] × [0, 30],  t ∈ (0, 20]"


# --- cells_per_axis -----------------------------------------------------------------------------------------------

def test_one_count_is_every_axis_and_a_map_is_one_per_axis():
    assert cells_per_axis(49, 2) == [49, 49]
    assert cells_per_axis(8, 3) == [8, 8, 8]
    assert cells_per_axis({"x": 80, "y": 40}, 2) == [80, 40]
    assert cells_per_axis({"y": 40, "x": 80}, 2) == [80, 40]            # in axis order, whatever the order written
    assert cells_per_axis({"x": 4, "y": 5, "z": 6}, 3) == [4, 5, 6]
    assert cells_per_axis({"x": 80, "y": 40}) == [80, 40]               # no dimension needed for a map


def test_a_key_reads_cells_per_axis_as_x_by_y():
    assert cells_text(49) == "49" and cells_text({"x": 80, "y": 40}) == "80x40" and cells_text({"y": 2, "x": 1}) == "1x2"


@pytest.mark.parametrize("value, dim, message", [
    ({"x": 80}, 2, "a 2D box has the axes x, y, and this gives x"),
    ({"x": 80, "y": 40, "z": 10}, 2, "a 2D box has the axes x, y, and this gives x, y, z"),
    ({"x": 80, "y": 40}, 3, "a 3D box has the axes x, y, z"),
    ({"x": 80, "w": 40}, 2, "no axis w; the axes are x, y, z"),
    ({"x": 80, "z": 40}, None, "none skipped"),
    ({"x": 0, "y": 40}, 2, r"mesh\.cells\.x: 0 is not a whole number of cells, at least 1"),
    ({"x": 80, "y": 4.5}, 2, r"mesh\.cells\.y: 4\.5 is not a whole number"),
    ({"x": 80, "y": "many"}, 2, "'many' is not a whole number"),
    ([80, 40], 2, "a list is a study, and has no axes"),
    (49, None, "needs the box's dimension"),
    (-3, 2, "not a whole number of cells, at least 1"),
])
def test_cells_that_do_not_make_sense_are_refused(value, dim, message):
    with pytest.raises(CellsError, match=message):
        cells_per_axis(value, dim)


@pytest.mark.parametrize("domain, dim", [
    (BOX, 2), ("Ω = [0, 1]^3", 3), ("Ω = [0, 1]^2,  t ∈ (0, 2]", 2), ("Ω = [0, 1] × [0, 2] × [0, 3]", 3),
    ("Ω = [0, 25] × [0, 1]", 2), ({"kind": "box", "lower": [0, 0], "upper": [1, 1]}, 2),
    ({"kind": "file", "format": "gmsh_geo"}, None), ("a disc", None),
])
def test_the_dimension_of_a_box_is_read_without_the_symbolic_front_end(domain, dim):
    assert domain_dim(domain) == dim


# --- the spec -----------------------------------------------------------------------------------------------------

@pytest.fixture(scope="module")
def dam():
    return load_ymf(DAM3BUMPS)


def test_dam3bumps_states_the_decks_80_by_40_cells(dam):
    # test/SWFlow/dam3Bumps.py, refinement 4: nnx = 5·2⁴ + 1 = 81, nny = (nnx − 1)/2 + 1 = 41; C5 could not say it
    (entry,) = dam["solution_paths"]["discretizations"]
    nnx = (6 - 1) * 2 ** 4 + 1
    nny = (nnx - 1) // 2 + 1
    assert cells_per_axis(entry["mesh"]["cells"], domain_dim(dam["model"]["strong_form"]["domain"])) == [nnx - 1, nny - 1]
    assert reached(dam)[0] == "discretized"


def test_it_is_one_realization_with_its_own_key(dam):
    (entry,) = dam["solution_paths"]["discretizations"]
    (realization,) = closure.realizations(entry, dam["model"]["strong_form"]["domain"])
    assert realization == {"cells": {"x": 80, "y": 40}, "levels": 1}
    (key,) = [k for k, _, _ in closure.planned_outputs(dam)]
    assert key.startswith("convex_limited/cells=80x40/")


def test_80_by_40_and_40_by_80_and_80_by_80_are_three_inputs(dam):
    seen = set()
    for cells in ({"x": 80, "y": 40}, {"x": 40, "y": 80}, {"x": 80, "y": 80}, 80):
        doc = copy.deepcopy(dam)
        doc["solution_paths"]["discretizations"][0]["mesh"]["cells"] = cells
        seen.add(closure.input_digest(doc, "convex_limited", {"cells": cells, "levels": 1}))
    assert len(seen) == 4        # even {x: 80, y: 80} is not 80: one says what it means, the other leaves it to the box


def test_a_list_of_maps_is_a_study_of_meshes(dam):
    doc = copy.deepcopy(dam)
    entry = doc["solution_paths"]["discretizations"][0]
    entry["mesh"]["cells"] = [{"x": 80, "y": 40}, {"x": 160, "y": 80}, 20]
    assert closure.realizations(entry, doc["model"]["strong_form"]["domain"]) == [
        {"cells": {"x": 80, "y": 40}, "levels": 1}, {"cells": {"x": 160, "y": 80}, "levels": 1}, {"cells": 20, "levels": 1}]
    check(doc)
    assert [k.split("/")[1] for k, _, _ in closure.planned_outputs(doc)] == ["cells=80x40", "cells=160x80", "cells=20"]


def test_an_older_spec_keeps_its_keys():
    # a spec with one count, or a study of counts, has the key it always had
    doc = load_ymf(SHALLOW_WATER / "lake_at_rest.ymf")
    assert [k.split("/")[1] for k, _, _ in closure.planned_outputs(doc)] == ["cells=50", "cells=100", "cells=200"]


def test_cells_and_h_are_still_alternatives(dam):
    doc = copy.deepcopy(dam)
    doc["solution_paths"]["discretizations"][0]["mesh"]["h"] = 2.0
    with pytest.raises(ValueError, match="mesh gives both cells and h"):
        closure.realizations(doc["solution_paths"]["discretizations"][0], doc["model"]["strong_form"]["domain"])


# --- checked on load ----------------------------------------------------------------------------------------------

@pytest.mark.parametrize("cells, message", [
    ({"x": 80}, "discretization convex_limited mesh.cells: a 2D box has the axes x, y"),
    ({"x": 80, "y": 40, "z": 5}, "a 2D box has the axes x, y"),
    ({"x": 80, "q": 40}, "no axis q"),
    ([{"x": 80, "y": 40}, {"x": 80}], "a 2D box has the axes x, y"),
])
def test_a_bad_cells_is_refused_on_load_and_names_the_discretization(dam, cells, message):
    doc = copy.deepcopy(dam)
    doc["solution_paths"]["discretizations"][0]["mesh"]["cells"] = cells
    with pytest.raises(VocabularyError, match=message):
        check(doc)


def test_the_strict_schema_refuses_a_cells_of_the_wrong_shape(tmp_path):
    for name in ("shallow_water.ymf", "proteus_numerics.ymf"):
        (tmp_path / name).write_text((SHALLOW_WATER / name).read_text(encoding="utf-8"), encoding="utf-8")
    text = DAM3BUMPS.read_text(encoding="utf-8").replace("cells: {x: 80, y: 40}", "cells: {x: eighty, y: 40}")
    (tmp_path / "bad.ymf").write_text(text, encoding="utf-8")
    with pytest.raises(Exception, match="eighty|integer"):
        load_ymf(tmp_path / "bad.ymf")


# --- the archive --------------------------------------------------------------------------------------------------

def test_cells_per_axis_survive_the_archive_and_the_key_is_the_same(tmp_path):
    realization = {"cells": {"x": 80, "y": 40}, "levels": 1}
    spec, again, digest = roundtrip(tmp_path, DAM3BUMPS, "convex_limited", realization)
    assert again["solution_paths"]["discretizations"][0]["mesh"]["cells"] == {"x": 80, "y": 40}
    assert [k for k, _, _ in closure.planned_outputs(again)] == [k for k, _, _ in closure.planned_outputs(spec)]
