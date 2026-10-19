"""ymf.discretizations: what each method needs before it can plan an output,
and ymf.state.branches: how far each discretization has got."""

from pathlib import Path

import pytest

from ymf import closure
from ymf.discretizations import METHODS, Missing, method_of, missing

POISSON = (Path(__file__).parent.parent / "examples" / "poisson.ymf").read_text(encoding="utf-8")
MESH_HINT = ("give mesh: {h: ...} (or {cells: ...} on a box), or --h (--cells) "
             "on the command line")


def test_mesh_is_the_default_method():
    assert method_of({"name": "P1"}) is METHODS["mesh"]
    assert method_of({"name": "P1", "method": "particle"}) is METHODS["particle"]
    with pytest.raises(ValueError, match="method 'fem' is not one of mesh, neural, particle"):
        method_of({"name": "P1", "method": "fem"})


def test_a_mesh_entry_needs_one_mesh_size():
    assert missing({"name": "P1"}) == [Missing("a mesh", MESH_HINT)]
    assert missing({"name": "P1", "mesh": {"levels": 2}}) == [Missing("a mesh", MESH_HINT)]
    assert missing({"name": "P1", "mesh": {"cells": 4, "h": 0.1}})[0].what == "one mesh size"
    assert missing({"name": "P1", "mesh": {"h": [0.1, 0.05]}}) == []


def test_neural_and_particle_name_each_missing_field():
    lacking = missing({"name": "N", "method": "neural", "neural": {"width": 32}})
    assert [m.what for m in lacking] == ["neural.depth", "neural.points"]
    assert lacking[0].hint == "give neural: {depth: ...}: the number of hidden layers"
    assert [m.what for m in missing({"name": "S", "method": "particle"})] == [
        "particle.spacing", "particle.kernel"]


def test_a_complete_entry_of_an_unimplemented_method_says_so():
    entry = {"name": "N", "method": "neural",
             "neural": {"width": 32, "depth": 4, "points": [1000, 4000]}}
    (only,) = missing(entry)
    assert only.what == "an implementation"
    assert "method neural is a schema only so far" in only.hint
    assert closure.realizations(entry) == []


def test_a_block_of_another_method_is_named():
    lacking = missing({"name": "N", "method": "neural", "mesh": {"cells": 4},
                       "neural": {"width": 32, "depth": 4, "points": 100}})
    assert lacking[0].hint == "mesh: belongs to method mesh, and this entry's method is neural"


def test_mesh_keys_read_as_they_always_have():
    mesh = METHODS["mesh"]
    assert mesh.key_parts({"cells": 4, "levels": 1}) == ["cells=4"]
    assert mesh.key_parts({"h": 0.1, "levels": 2, "dt": 0.5}) == ["h=0.1", "levels=2", "dt=0.5"]
    assert METHODS["neural"].key_parts({"width": 32, "depth": 4, "points": 100}) == [
        "width=32", "depth=4", "points=100"]


# --- through a document --------------------------------------------------------------

pytest.importorskip("strictyaml")

from ymf.archive import add_uniform_step, data_item, geometry, new_domain, topology  # noqa: E402
from ymf.compose import load_composed  # noqa: E402
from ymf.state import branches, reached  # noqa: E402


def _output(realization):
    domain = new_domain()
    nodes = data_item([2, 2], values=[0.0, 0.0, 1.0, 1.0])
    vertices = data_item([2], values=[0, 1], data_type="Int", precision=4)
    add_uniform_step(domain, 0.0, topology("Polyvertex", 2, vertices, nodes_per_element=1),
                     geometry(nodes, geometry_type="XY"))
    return {"realization": realization, "approximation": domain}


def _without_the_second_mesh(text=POISSON):
    """examples/poisson.ymf with no mesh on P2_quadratic."""
    head, tail = text.split('    - name: "P2_quadratic"\n', 1)
    tail = "".join(l for l in tail.splitlines(True) if not l.lstrip().startswith("mesh:"))
    return head + '    - name: "P2_quadratic"\n' + tail


def test_each_branch_reaches_its_own_state(tmp_path):
    path = tmp_path / "poisson.ymf"
    path.write_text(_without_the_second_mesh(), encoding="utf-8")
    doc = load_composed(path)
    doc["outputs"] = {k: _output(r) for k, _, r in closure.planned_outputs(doc)}
    p1, p2 = branches(doc)
    assert (p1.name, p1.state, p1.missing) == ("P1_linear", "closed", ())
    assert len(p1.planned) == 3
    assert (p2.name, p2.state) == ("P2_quadratic", "well-posed")
    assert p2.reasons == ("P2_quadratic is missing a mesh: " + MESH_HINT,)
    # the document is as far as its least advanced branch, and says which
    assert reached(doc) == ("well-posed", p2.reasons[0])


def test_a_spec_may_name_a_method_ymf_cannot_plan_yet(tmp_path):
    text = POISSON.replace(
        '    - name: "P2_quadratic"\n',
        '    - name: "P2_quadratic"\n'
        '      method: particle\n'
        '      particle: {spacing: [0.1, 0.05], kernel: cubic_spline}\n', 1)
    text = _without_the_second_mesh(text)
    path = tmp_path / "poisson.ymf"
    path.write_text(text, encoding="utf-8")
    doc = load_composed(path)                     # validates: well-posed
    state, why = reached(doc)
    assert state == "well-posed"
    assert why.startswith("P2_quadratic is missing an implementation: method particle")
