"""ymf.state: how far a document has got, declared and checked."""

import re
from pathlib import Path

import pytest

pytest.importorskip("strictyaml")

from ymf import closure  # noqa: E402
from ymf.archive import add_uniform_step, data_item, geometry, new_domain, topology  # noqa: E402
from ymf.compose import YmfCompositionError, load_composed  # noqa: E402
from ymf.state import STATES, YmfStateError, reached  # noqa: E402

POISSON = (Path(__file__).parent.parent / "examples" / "poisson.ymf").read_text(encoding="utf-8")


def _output(realization):
    domain = new_domain()
    nodes = data_item([2, 2], values=[0.0, 0.0, 1.0, 1.0])
    vertices = data_item([2], values=[0, 1], data_type="Int", precision=4)
    add_uniform_step(domain, 0.0, topology("Polyvertex", 2, vertices, nodes_per_element=1),
                     geometry(nodes, geometry_type="XY"))
    return {"realization": realization, "approximation": domain}


def _write(tmp_path, text, name="poisson.ymf"):
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


def _declaring(state, text=POISSON):
    return "ymf_state: %s\n%s" % (state, text)


def _without_meshes(text=POISSON):
    return "".join(line for line in text.splitlines(True)
                   if not line.lstrip().startswith("mesh:"))


def test_the_states_in_order():
    assert STATES == ("abstract", "well-posed", "discretized", "closed")


def test_each_state_is_reached_when_its_content_is(tmp_path):
    meshed = load_composed(_write(tmp_path, POISSON))
    assert reached(meshed) == ("discretized",
                               "6 of 6 planned outputs are missing: " +
                               ", ".join(k for k, _, _ in closure.planned_outputs(meshed)))
    assert reached(meshed, well_posed=False)[0] == "abstract"
    unmeshed = load_composed(_write(tmp_path, _without_meshes(), "bare.ymf"))
    state, why = reached(unmeshed)
    assert state == "well-posed"
    assert why == "; ".join("%s is missing a mesh: give mesh: {h: ...} (or {cells: ...} on "
                            "a box), or --h (--cells) on the command line" % n
                            for n in ("P1_linear", "P2_quadratic"))


def test_a_declared_state_the_content_reaches_loads(tmp_path):
    for state in ("well-posed", "discretized"):
        doc = load_composed(_write(tmp_path, _declaring(state)))
        assert doc["ymf_state"] == state


def test_a_declared_state_beyond_the_content_is_refused(tmp_path):
    with pytest.raises(YmfStateError, match="declares ymf_state: closed, but it is only "
                                             "discretized: 6 of 6 planned outputs"):
        load_composed(_write(tmp_path, _declaring("closed")))
    with pytest.raises(YmfStateError, match="only well-posed: P1_linear is missing a mesh"):
        load_composed(_write(tmp_path, _declaring("discretized", _without_meshes())))


def test_an_archive_records_the_state_it_reaches(tmp_path):
    spec, _, _ = closure.load(_write(tmp_path, POISSON))
    outputs = {k: _output(r) for k, _, r in closure.planned_outputs(spec)}
    archive = tmp_path / "poisson.archive.ymf"
    closure.write(archive, spec, outputs)
    assert "\nymf_state: closed\n" in archive.read_text(encoding="utf-8")
    closure.load(archive)                                  # and it holds

    partial = tmp_path / "partial.archive.ymf"
    closure.write(partial, spec, dict(list(outputs.items())[:2]))
    assert "\nymf_state: discretized\n" in partial.read_text(encoding="utf-8")


def test_an_archive_that_no_longer_reaches_its_state(tmp_path):
    spec, _, _ = closure.load(_write(tmp_path, POISSON))
    archive = tmp_path / "poisson.archive.ymf"
    closure.write(archive, spec, {k: _output(r) for k, _, r in closure.planned_outputs(spec)})
    # edit the input inside the archive: every output is keyed to the old one
    text = archive.read_text(encoding="utf-8")
    edited = re.sub(r"name: Poisson unit square", "name: Poisson edited", text, count=1)
    assert edited != text
    archive.write_text(edited, encoding="utf-8")
    with pytest.raises(YmfStateError, match="declares ymf_state: closed, but it is only "
                                             "discretized"):
        closure.load(archive)
    # a runner bringing it up to date loads it anyway, and recomputes
    spec, kept, dropped = closure.load(archive, check_state=False)
    assert kept == {} and len(dropped) == 6


def test_an_abstract_file_stays_abstract_through_composition(tmp_path):
    _write(tmp_path, _declaring("abstract", _without_meshes()), "base.ymf")
    child = _write(tmp_path, "ymf_state: abstract\nextends: base.ymf\n", "child.ymf")
    assert load_composed(child)["ymf_state"] == "abstract"
    # the state describes the file: a child that declares none is checked
    # as a full document, whatever its parent declared
    plain = _write(tmp_path, "extends: base.ymf\n", "plain.ymf")
    doc = load_composed(plain)
    assert "ymf_state" not in doc


def test_the_old_file_level_kind_says_what_replaced_it(tmp_path):
    with pytest.raises(YmfCompositionError, match="is now 'ymf_state': write "
                                                   "'ymf_state: abstract' for 'kind: model'"):
        load_composed(_write(tmp_path, "kind: model\n" + _without_meshes()))
