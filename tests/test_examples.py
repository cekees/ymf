"""The examples are documentation, so they are run here to keep them true."""

import subprocess
import sys
from pathlib import Path

import pytest

from ymf.schema import load_ymf

EXAMPLES_DIR = Path(__file__).parent.parent / "examples"
SPECS = sorted(EXAMPLES_DIR.glob("*.yaml"))


def run(script, *args):
    return subprocess.run([sys.executable, str(EXAMPLES_DIR / script), *map(str, args)],
                          capture_output=True, text=True, check=True).stdout


@pytest.mark.parametrize("spec", SPECS, ids=[p.name for p in SPECS])
def test_every_example_spec_validates(spec):
    assert load_ymf(spec)["Problem"]["name"]


def test_validate_spec_walks_every_example_and_shows_the_unit_error():
    out = run("validate_spec.py")
    for spec in SPECS:
        assert spec.name in out
    assert "schema: INVALID" not in out
    assert "dimensionless numbers: {'Re': 40.0}" in out
    assert "incompatible with 'T''s declared units 'K'" in out


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
