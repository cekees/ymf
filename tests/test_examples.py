"""The examples are documentation, so they are run here to keep them true."""

import os
import subprocess
import sys
from pathlib import Path

import pytest

from ymf.schema import load_ymf

EXAMPLES_DIR = Path(__file__).parent.parent / "examples"
sys.path.insert(0, str(EXAMPLES_DIR))
from validate_spec import standalone  # noqa: E402

SPECS = [p for p in sorted(EXAMPLES_DIR.rglob("*.ymf")) if standalone(p)]


def run(script, *args, env=None):
    return subprocess.run([sys.executable, str(EXAMPLES_DIR / script), *map(str, args)],
                          capture_output=True, text=True, check=True, env=env).stdout


@pytest.mark.parametrize("spec", SPECS, ids=[p.name for p in SPECS])
def test_every_example_spec_validates(spec):
    doc = load_ymf(spec)
    # a model is not a problem yet, so only problems carry a name
    assert doc.get("kind") == "model" or doc["Problem"]["name"]


def test_validate_spec_walks_every_example_and_shows_the_unit_error():
    out = run("validate_spec.py")
    for spec in SPECS:
        assert spec.name in out
    assert "composed from: navier_stokes_model.ymf <- plane_poiseuille.ymf <- plane_poiseuille_re100.ymf" in out
    assert "override Problem.strong_form.coefficients.μ.value: '1.0' -> '0.01' (plane_poiseuille_re100.ymf)" in out
    assert "schema: INVALID" not in out
    assert "dimensionless numbers: {'Re': 40.0}" in out
    assert "incompatible with 'T''s declared units 'K'" in out


def test_the_refined_example_inherits_every_output_its_change_leaves_alone(tmp_path):
    import shutil
    from ymf import closure
    here = tmp_path / "navier_stokes"
    shutil.copytree(EXAMPLES_DIR / "navier_stokes", here,
                    ignore=shutil.ignore_patterns("*.archive.ymf", "*.h5", "*.xmf"))
    spec, _, _ = closure.load(here / "plane_poiseuille.ymf")
    planned = list(closure.planned_outputs(spec))
    assert len(planned) == 8
    closure.write(here / "plane_poiseuille.archive.ymf", spec,
                  {key: {"realization": r} for key, _, r in planned})
    refined, kept, dropped = closure.load(here / "plane_poiseuille_refined.ymf")
    assert len(kept) == 8 and dropped == {}
    new = [k for k, _, _ in closure.planned_outputs(refined) if k not in kept]
    assert len(new) == 1 and new[0].startswith("equal_order_p1/cells=128/")


@pytest.fixture(scope="module")
def archive_dir(tmp_path_factory):
    pytest.importorskip("numpy")
    pytest.importorskip("h5py")
    outdir = tmp_path_factory.mktemp("archives")
    run("write_archive.py", "--outdir", outdir, "-n", "8")
    return outdir


@pytest.mark.parametrize("name", ["heat", "heat_split"])
def test_write_archive_produces_ymf_h5_and_xmf(archive_dir, name):
    for suffix in (".ymf", ".h5", ".xmf"):
        assert (archive_dir / (name + suffix)).is_file()


@pytest.mark.parametrize("name", ["heat", "heat_split"])
def test_read_archive_recovers_kappa_and_matches_the_xmf(archive_dir, name):
    out = run("read_archive.py", archive_dir / (name + ".ymf"))
    assert "kappa fitted from the archived fields: 0.001\n" in out
    assert "holds the same domain and problem: True" in out


@pytest.fixture(scope="module")
def inline_dir(tmp_path_factory):
    """Archives written with --inline, with h5py made unimportable."""
    pytest.importorskip("numpy")
    blocker = tmp_path_factory.mktemp("block")
    (blocker / "h5py.py").write_text('raise ImportError("h5py is blocked for this test")\n')
    env = dict(os.environ, PYTHONPATH=os.pathsep.join(
        [str(blocker)] + [p for p in os.environ.get("PYTHONPATH", "").split(os.pathsep) if p]))
    outdir = tmp_path_factory.mktemp("inline")
    run("write_archive.py", "--outdir", outdir, "-n", "8", "--inline", env=env)
    return outdir, env


@pytest.mark.parametrize("name", ["heat", "heat_split"])
def test_inline_archives_are_self_contained_and_need_no_h5py(inline_dir, name):
    outdir, env = inline_dir
    assert sorted(p.name for p in outdir.iterdir()) == [
        "heat.xmf", "heat.ymf", "heat_split.xmf", "heat_split.ymf"]
    out = run("read_archive.py", outdir / (name + ".ymf"), env=env)
    assert "kappa fitted from the archived fields: 0.001\n" in out
    assert "holds the same domain and problem: True" in out
