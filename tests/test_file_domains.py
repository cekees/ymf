"""`file` domains (C27 Phase 5, D2 and D3).

A file domain is carried whole for the solver: a mesh file (2dm, 3dm, hex,
tetgen, ymf) fixes the mesh, so a discretization gives only its `levels` of
uniform refinement; a geometry (gmsh_geo) is meshed by `h`. Its `bounds` are
optional: without them the solver takes them from the mesh. Boundary
regions select its boundary tags; periodic conditions are not carried.
"""

import copy
from pathlib import Path

import pytest

pytest.importorskip("strictyaml")

from ymf.domain import YmfDomainError, validate_domain_block  # noqa: E402
from ymf.references import YmfReferenceError, resolve  # noqa: E402
from ymf.schema import load_ymf  # noqa: E402
from ymf.state import reached  # noqa: E402

EXAMPLES = Path(__file__).parent.parent / "examples"

CUBE = {"kind": "file", "format": "tetgen", "path": "cube", "dimension": 3,
        "boundary_tags": {"bottom": 1, "top": 2, "front": 3, "back": 4, "left": 5,
                          "right": 6}}


def _on_file(domain, mesh):
    doc = load_ymf(EXAMPLES / "poisson.ymf")
    strong = doc["model"]["strong_form"]
    strong["domain"] = copy.deepcopy(domain)
    strong["boundary_regions"] = [{"name": "fixed", "tags": ["bottom", "top"]},
                                  {"name": "sides", "tags": [3, 4, 5, 6]}]
    strong["boundary_conditions"] = [
        {"region": "fixed", "variable": "u", "type": "dirichlet", "value": 0.0},
        {"region": "sides", "variable": "u", "type": "neumann", "value": 0.0}]
    doc["solution_paths"]["analytical"] = []
    doc.pop("vvuq", None)
    for d in doc["solution_paths"]["discretizations"]:
        d["mesh"] = copy.deepcopy(mesh)
    return doc


def test_a_mesh_file_plans_its_levels():
    from ymf.closure import planned_outputs
    doc = _on_file(CUBE, {"levels": [1, 2]})
    assert reached(doc)[0] == "discretized"
    keys = sorted(k for k, _, _ in planned_outputs(doc))
    assert [k.rsplit("/", 1)[0] for k in keys] == [
        "P1_linear/levels=1", "P1_linear/levels=2", "P2_quadratic/levels=1",
        "P2_quadratic/levels=2"]


def test_a_mesh_file_takes_no_mesh_size():
    state, why = reached(_on_file(CUBE, {"h": 0.1}))
    assert state == "well-posed"
    assert "already a mesh; give mesh: {levels: ...} only" in why


def test_a_geometry_file_is_meshed_by_h():
    from ymf.closure import planned_outputs
    geo = dict(CUBE, format="gmsh_geo", path="cube.geo")
    keys = [k for k, _, _ in planned_outputs(_on_file(geo, {"h": [0.2, 0.1]}))]
    assert keys[0].startswith("P1_linear/h=0.2/")


def test_the_tags_resolve_through_the_files_boundary_tags():
    pytest.importorskip("sympy")
    from ymf.symbolic import adr_problem
    problem = adr_problem(_on_file(CUBE, {"levels": 1}))
    assert problem["domain"] == {"lower": None, "upper": None,
                                 "region": problem["domain"]["region"]}
    assert problem["dirichlet"][0]["flags"] == [1, 2]
    assert problem["neumann"][0]["flags"] == [3, 4, 5, 6]


def test_bounds_set_the_box():
    pytest.importorskip("sympy")
    from ymf.symbolic import adr_problem
    cube = dict(CUBE, bounds={"lower": [0, 0, 0], "upper": [1, 2, 3]})
    domain = adr_problem(_on_file(cube, {"levels": 1}))["domain"]
    assert domain["lower"] == [0.0, 0.0, 0.0] and domain["upper"] == [1.0, 2.0, 3.0]


@pytest.mark.parametrize("bounds, message", [
    ({"lower": [0, 0], "upper": [1, 1, 1]}, "same, nonzero, number"),
    ({"lower": [0, 0, 1], "upper": [1, 1, 1]}, r"upper\[2\] = 1 is not above"),
    ({"lower": [0, 0], "upper": [1, 1]}, "dimension 3 but bounds of 2"),
])
def test_bounds_are_checked(bounds, message):
    with pytest.raises(YmfDomainError, match=message):
        validate_domain_block(dict(CUBE, bounds=bounds))


def test_a_file_domain_is_not_periodic():
    pytest.importorskip("sympy")
    from ymf.symbolic import SymbolicError, adr_problem
    doc = _on_file(dict(CUBE, bounds={"lower": [0, 0, 0], "upper": [1, 1, 1]}), {"levels": 1})
    strong = doc["model"]["strong_form"]
    strong["boundary_regions"].append({"name": "ends", "where": "x = 0 or x = 1"})
    strong["boundary_conditions"].append({"region": "ends", "variable": "u",
                                          "type": "periodic"})
    resolve(doc)
    with pytest.raises(SymbolicError, match="a file domain cannot be periodic"):
        adr_problem(doc)


def test_a_file_domain_loads_with_bounds(tmp_path):
    text = (EXAMPLES / "poisson.ymf").read_text(encoding="utf-8").replace(
        '''    domain: "Ω = [0, 1] × [0, 1]"''', '''    domain:
      kind: file
      format: tetgen
      path: cube
      bounds: {lower: [0, 0, 0], upper: [1, 1, 1]}
      boundary_tags: {walls: 1}''').replace('''        where: "∂Ω"''', '''        tags: [walls]''').replace(
        "mesh: {cells: [4, 8, 16]}", "mesh: {levels: [1, 2]}")
    (tmp_path / "p.ymf").write_text(text, encoding="utf-8")
    doc = load_ymf(tmp_path / "p.ymf")
    assert doc["model"]["strong_form"]["domain"]["bounds"]["upper"] == [1.0, 1.0, 1.0]
    assert reached(doc)[0] == "discretized"
