"""`discrete` domains: a mesh as ymf stores one, solved as it stands (C27
Phase 5, D4; C7.1: a mesh is a type of domain).

Its bounds are its nodes' when the document holds them, and the solver's to
take when they are in a file; its boundary regions select the tags of its
facets (#6); a discretization gives only its levels; it is not periodic.
The solver reads it through Proteus's MeshYmfDomain (cekees/proteus#186).
"""

import copy
from pathlib import Path

import pytest

pytest.importorskip("strictyaml")
pytest.importorskip("sympy")

from ymf.domain import validate_domain_block  # noqa: E402
from ymf.schema import load_ymf  # noqa: E402
from ymf.state import reached  # noqa: E402
from ymf.symbolic import SymbolicError, adr_problem  # noqa: E402

EXAMPLES = Path(__file__).parent.parent / "examples"


def _item(values, dims, data_type="Float"):
    return {"data_type": data_type, "precision": 8 if data_type == "Float" else 4,
            "dimensions": dims, "values": values}


#: the square [0, 2] x [0, 1] as two triangles, its four edges tagged
SQUARE = {
    "kind": "discrete",
    "topology": {"type": "Triangle", "number_of_elements": 2,
                 "data_item": _item([0, 1, 2, 0, 2, 3], [2, 3], "Int")},
    "geometry": {"type": "XY", "data_item": _item([0, 0, 2, 0, 2, 1, 0, 1], [4, 2])},
    "facets": {"topology": {"type": "Polyline", "number_of_elements": 4,
                            "data_item": _item([0, 1, 1, 2, 2, 3, 3, 0], [4, 2], "Int")},
               "tags": {"data_item": _item([1, 2, 3, 4], [4], "Int")},
               "view": "global"},
    "boundary_tags": {"bottom": 1, "right": 2, "top": 3, "left": 4},
}


def _on(domain, mesh):
    doc = load_ymf(EXAMPLES / "poisson.ymf")
    strong = doc["model"]["strong_form"]
    strong["domain"] = copy.deepcopy(domain)
    strong["boundary_regions"] = [{"name": "fixed", "tags": ["bottom", "left"]},
                                  {"name": "rest", "tags": [2, 3]}]
    strong["boundary_conditions"] = [
        {"region": "fixed", "variable": "u", "type": "dirichlet", "value": 0.0},
        {"region": "rest", "variable": "u", "type": "neumann", "value": 0.0}]
    doc["solution_paths"]["analytical"] = []
    doc.pop("vvuq", None)
    for d in doc["solution_paths"]["discretizations"]:
        d["mesh"] = copy.deepcopy(mesh)
    return doc


def test_a_mesh_is_a_domain():
    validate_domain_block(SQUARE)


def test_its_bounds_are_its_nodes():
    problem = adr_problem(_on(SQUARE, {"levels": 1}))
    assert problem["domain"]["lower"] == [0.0, 0.0]
    assert problem["domain"]["upper"] == [2.0, 1.0]
    assert problem["domain"]["region"]["kind"] == "discrete"
    assert problem["dirichlet"][0]["flags"] == [1, 4]
    assert problem["neumann"][0]["flags"] == [2, 3]


def test_nodes_in_a_file_leave_the_bounds_to_the_solver():
    square = copy.deepcopy(SQUARE)
    square["geometry"]["data_item"] = {"format": "HDF", "data_type": "Float",
                                       "precision": 8, "dimensions": [4, 2],
                                       "data": "square.h5:/nodes"}
    problem = adr_problem(_on(square, {"levels": 1}))
    assert problem["domain"]["lower"] is None
    assert problem["adr"]["dim"] == 2                     # from the geometry's type
    square["bounds"] = {"lower": [0, 0], "upper": [2, 1]}
    assert adr_problem(_on(square, {"levels": 1}))["domain"]["upper"] == [2.0, 1.0]


def test_a_discretization_refines_it_by_levels():
    from ymf.closure import planned_outputs
    doc = _on(SQUARE, {"levels": [1, 2]})
    assert reached(doc)[0] == "discretized"
    assert sorted(k.rsplit("/", 1)[0] for k, _, _ in planned_outputs(doc))[:2] == [
        "P1_linear/levels=1", "P1_linear/levels=2"]
    assert "already a mesh" in reached(_on(SQUARE, {"cells": 4}))[1]


def test_a_discrete_domain_is_not_periodic():
    doc = _on(SQUARE, {"levels": 1})
    strong = doc["model"]["strong_form"]
    strong["boundary_regions"].append({"name": "ends", "where": "x = 0 or x = 2"})
    strong["boundary_conditions"].append({"region": "ends", "variable": "u",
                                          "type": "periodic"})
    with pytest.raises(SymbolicError, match="a discrete domain cannot be periodic"):
        adr_problem(doc)
