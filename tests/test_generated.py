"""ymf.generated: files a tool wrote, and whether they still hold."""

import copy
from pathlib import Path

import pytest

pytest.importorskip("strictyaml")

from ymf import closure, generated
from ymf.archive import read_document
from ymf.compose import load_composed

EXAMPLES = Path(__file__).parent.parent / "examples"
BY = "pro weak (test)"


def _weak_file(tmp_path):
    """examples/poisson.ymf, plus a generated file adding one weak form."""
    spec = tmp_path / "poisson.ymf"
    spec.write_text((EXAMPLES / "poisson.ymf").read_text(encoding="utf-8"), encoding="utf-8")
    parent = load_composed(spec)
    weak = copy.deepcopy(parent["model"]["weak_forms"][0])
    weak.update(label="galerkin_generated", provenance="proteus_derived")
    path = tmp_path / "poisson.weak.ymf"
    record = generated.write(path, "poisson.ymf", BY, "model_without_weak_forms",
                             parent, {"model": {"weak_forms": [weak]}})
    return spec, path, record


def test_a_written_file_is_current_and_composes(tmp_path):
    _, path, record = _weak_file(tmp_path)
    assert record["by"] == BY and len(record["from"]) == len(record["content"]) == 12
    status = generated.check(path)
    assert status.current and status.reasons == []
    doc = load_composed(path)
    assert [w["label"] for w in doc["model"]["weak_forms"]] == ["global", "galerkin_generated"]
    assert doc["generated"] == record                  # the top file's, carried


def test_an_upstream_edit_it_read_makes_it_stale(tmp_path):
    spec, path, _ = _weak_file(tmp_path)
    text = spec.read_text(encoding="utf-8")
    spec.write_text(text.replace('name: "Poisson unit square"', 'name: "Poisson, edited"'),
                    encoding="utf-8")
    status = generated.check(path)
    assert status.stale and not status.edited and not status.current
    assert "poisson.ymf changed in what pro weak (test) read" in status.reasons[0]


def test_an_upstream_edit_it_did_not_read_leaves_it_current(tmp_path):
    spec, path, _ = _weak_file(tmp_path)
    text = spec.read_text(encoding="utf-8")
    old = 'a(u,v) = ∫Ω ∇u·∇v dx'
    assert old in text
    spec.write_text(text.replace(old, 'a(u,v) = ∫Ω ∇v·∇u dx'), encoding="utf-8")
    assert generated.check(path).current           # pro weak does not read weak forms


def test_an_edit_to_what_it_wrote_is_reported(tmp_path):
    _, path, _ = _weak_file(tmp_path)
    text = path.read_text(encoding="utf-8")
    path.write_text(text.replace("proteus_derived", "human_edited"), encoding="utf-8")
    status = generated.check(path)
    assert status.edited and not status.stale
    assert "content was edited since pro weak (test) wrote it" in status.reasons[0]


def test_a_removal_counts_as_content(tmp_path):
    _, path, _ = _weak_file(tmp_path)
    with path.open("a", encoding="utf-8") as f:
        f.write("remove: ['model.physical_model.assumptions[steady]']\n")
    assert generated.check(path).edited


def test_it_is_not_inherited_and_archives_do_not_carry_it(tmp_path):
    _, path, _ = _weak_file(tmp_path)
    child = tmp_path / "poisson.disc.ymf"
    child.write_text("extends: poisson.weak.ymf\nmodel:\n  name: \"Poisson, child\"\n",
                     encoding="utf-8")
    doc = load_composed(child)
    assert "generated" not in doc
    assert generated.check(child).record is None
    archive = tmp_path / "poisson.weak.archive.ymf"
    closure.write(archive, load_composed(path), {})
    assert "generated" not in read_document(archive)


def test_what_a_step_reads_is_named_and_checked(tmp_path):
    with pytest.raises(ValueError, match="is not one of"):
        generated.slice_of({}, "everything")
    lone = tmp_path / "lone.ymf"
    lone.write_text("generated: {by: x, reads: model_and_analytical, from: a, content: b}\n"
                    "model:\n  name: \"x\"\n", encoding="utf-8")
    with pytest.raises(ValueError, match="extends nothing"):
        generated.check(lone)


def test_a_step_can_declare_the_state_it_reaches(tmp_path):
    from ymf.state import YmfStateError
    spec = tmp_path / "poisson.ymf"
    spec.write_text((EXAMPLES / "poisson.ymf").read_text(encoding="utf-8"), encoding="utf-8")
    parent = load_composed(spec)
    weak = copy.deepcopy(parent["model"]["weak_forms"][0])
    weak.update(label="galerkin_generated", provenance="proteus_derived")
    content = {"model": {"weak_forms": [weak]}}
    plain = generated.write(tmp_path / "plain.ymf", "poisson.ymf", BY,
                            "model_without_weak_forms", parent, content)
    path = tmp_path / "declared.ymf"
    record = generated.write(path, "poisson.ymf", BY, "model_without_weak_forms",
                             parent, content, state="discretized")
    assert record["content"] == plain["content"]      # the state is not content
    assert path.read_text(encoding="utf-8").splitlines()[1] == "ymf_state: discretized"
    assert generated.check(path).current
    assert load_composed(path)["ymf_state"] == "discretized"   # and the loader checked it
    generated.write(path, "poisson.ymf", BY, "model_without_weak_forms", parent, content,
                    state="closed")
    with pytest.raises(YmfStateError):
        load_composed(path)                            # declared, but not reached
    with pytest.raises(ValueError, match="not one of"):
        generated.write(path, "poisson.ymf", BY, "model_without_weak_forms", parent,
                        content, state="finished")
