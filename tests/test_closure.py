"""ymf.closure: realizations, output keys, and archives as spec + outputs."""

from pathlib import Path

import pytest

from ymf import closure
from ymf.archive import (YmfArchiveError, add_uniform_step, data_item, geometry,
                         new_domain, read_document, topology)

EXAMPLES = Path(__file__).parent.parent / "examples"
MESH = "      mesh: {cells: [4, 8], levels: 2}\n"


def spec_with_mesh(tmp_path, mesh=MESH, name="poisson.ymf"):
    """examples/poisson.ymf with a mesh study on P1_linear."""
    text = (EXAMPLES / "poisson.ymf").read_text(encoding="utf-8")
    anchor = '    - name: "P1_linear"\n'
    assert anchor in text
    text = text.replace(anchor, anchor + mesh, 1)
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


def fake_output(realization):
    """An output holding one inline, two-point grid."""
    domain = new_domain()
    nodes = data_item([2, 2], values=[0.0, 0.0, 1.0, 1.0])
    vertices = data_item([2], values=[0, 1], data_type="Int", precision=4)
    add_uniform_step(domain, 0.0, topology("Polyvertex", 2, vertices, nodes_per_element=1),
                     geometry(nodes, geometry_type="XY"))
    return {"realization": realization, "approximation": domain}


def test_realizations_are_the_cartesian_product_of_studied_values():
    entry = {"name": "P1", "mesh": {"cells": [4, 8], "levels": [1, 2]},
             "discretization": {"time_integration": "backward_euler", "dt": [0.1, 0.05]}}
    out = closure.realizations(entry)
    assert len(out) == 8
    assert out[0] == {"cells": 4, "levels": 1, "dt": 0.1}
    assert out[-1] == {"cells": 8, "levels": 2, "dt": 0.05}


def test_a_scalar_mesh_is_one_realization_and_no_mesh_is_none():
    assert closure.realizations({"mesh": {"cells": 8}}) == [{"cells": 8, "levels": 1}]
    assert closure.realizations({"name": "P1"}) == []


def test_realize_replaces_each_list_with_its_value():
    entry = {"name": "P1", "mesh": {"cells": [4, 8], "levels": 2}}
    assert closure.realize(entry, {"cells": 8, "levels": 2})["mesh"] == {"cells": 8, "levels": 2}
    assert entry["mesh"]["cells"] == [4, 8]           # untouched


def test_keys_name_the_realization_and_hash_the_slice(tmp_path):
    spec, _, _ = closure.load(spec_with_mesh(tmp_path))
    keys = [k for k, _, _ in closure.planned_outputs(spec)]
    assert [k.rsplit("/", 1)[0] for k in keys] == ["P1_linear/cells=4/levels=2",
                                                    "P1_linear/cells=8/levels=2"]
    assert all(len(k.rsplit("/", 1)[1]) == closure.KEY_DIGITS for k in keys)
    assert keys[0].rsplit("/", 1)[1] != keys[1].rsplit("/", 1)[1]


def test_a_key_ignores_other_discretizations_but_not_the_problem(tmp_path):
    spec, _, _ = closure.load(spec_with_mesh(tmp_path))
    r = {"cells": 4, "levels": 2}
    key = closure.output_key(spec, "P1_linear", r)
    spec["solution_paths"]["discretizations"][1]["solver"]["tolerance"] = 1e-6
    assert closure.output_key(spec, "P1_linear", r) == key
    spec["Problem"]["name"] = "renamed"                # prose counts too
    assert closure.output_key(spec, "P1_linear", r) != key


def test_an_archive_round_trips_and_its_embedded_spec_gives_the_same_keys(tmp_path):
    spec, _, _ = closure.load(spec_with_mesh(tmp_path))
    outputs = {k: fake_output(r) for k, _, r in closure.planned_outputs(spec)}
    archive = tmp_path / "poisson.archive.ymf"
    closure.write(archive, spec, outputs)
    document = read_document(archive)
    assert list(document)[:2] == ["ymf", "Problem"] and list(document)[-1] == "outputs"
    again, kept, dropped = closure.load(archive)
    assert set(kept) == set(outputs) and dropped == {}
    assert [k for k, _, _ in closure.planned_outputs(again)] == list(outputs)


def test_extending_an_archive_keeps_only_outputs_whose_input_is_unchanged(tmp_path):
    spec, _, _ = closure.load(spec_with_mesh(tmp_path))
    outputs = {k: fake_output(r) for k, _, r in closure.planned_outputs(spec)}
    closure.write(tmp_path / "poisson.archive.ymf", spec, outputs)
    child = tmp_path / "finer.ymf"
    # P2 gains a mesh; P1's slice is untouched, so its outputs carry over
    child.write_text(
        "extends: poisson.archive.ymf\n"
        "solution_paths:\n"
        "  discretizations:\n"
        "    - name: P2_quadratic\n"
        "      mesh: {cells: 4}\n", encoding="utf-8")
    merged, kept, dropped = closure.load(child)
    assert set(kept) == set(outputs) and dropped == {}
    assert merged["composition"]["sources"] == ["poisson.archive.ymf", "finer.ymf"]
    planned = [k for k, _, _ in closure.planned_outputs(merged)]
    assert set(outputs) < set(planned) and len(planned) == 3

    # changing the problem invalidates every inherited output
    child.write_text("extends: poisson.archive.ymf\nProblem:\n  name: Changed\n",
                     encoding="utf-8")
    _, kept, dropped = closure.load(child)
    assert kept == {} and set(dropped) == set(outputs)
    assert all("input changed" in why for why in dropped.values())


def test_a_kept_differing_rerun_is_matched_by_its_base_key(tmp_path):
    spec, _, _ = closure.load(spec_with_mesh(tmp_path, mesh="      mesh: {cells: 4}\n"))
    (key, _, r), = closure.planned_outputs(spec)
    closure.write(tmp_path / "a.archive.ymf", spec,
                  {key: fake_output(r), key + "~2": fake_output(r)})
    _, kept, _ = closure.load(tmp_path / "a.archive.ymf")
    assert set(kept) == {key, key + "~2"}


def test_a_document_from_a_newer_ymf_is_refused(tmp_path):
    path = tmp_path / "future.ymf"
    path.write_text("ymf: 999.0.0\nProblem: {name: x}\n", encoding="utf-8")
    with pytest.raises(YmfArchiveError, match="newer than this ymf"):
        read_document(path)


def test_ymf2xmf_writes_one_xmf_per_output(tmp_path):
    from ymf.cli import ymf2xmf_main
    from ymf.xdmf import read_xdmf
    spec, _, _ = closure.load(spec_with_mesh(tmp_path))
    outputs = {k: fake_output(r) for k, _, r in closure.planned_outputs(spec)}
    archive = tmp_path / "poisson.archive.ymf"
    closure.write(archive, spec, outputs)
    assert ymf2xmf_main([str(archive)]) == 0
    written = sorted(p.name for p in tmp_path.glob("*.xmf"))
    assert written == sorted("poisson.%s.xmf" % k.rsplit("/", 1)[1] for k in outputs)
    _, extra = read_xdmf(tmp_path / written[0])
    assert extra["output"] in outputs and "approximation" not in extra
    assert extra["Problem"]["name"] == spec["Problem"]["name"]
