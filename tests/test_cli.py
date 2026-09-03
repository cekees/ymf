"""Tests for the ymf2xdmf converter and its command line."""

import subprocess
import sys
from xml.etree.ElementTree import parse

import pytest

from ymf.archive import (
    YmfArchiveError,
    add_spatial_step,
    add_uniform_step,
    attribute,
    data_item,
    geometry,
    grid,
    new_domain,
    read_ymf,
    topology,
    write_ymf,
)
from ymf.cli import ymf2xdmf, ymf2xdmf_main


def make_domain(spatial=False, n_ranks=2):
    domain = new_domain("Mesh Spatial_Domain")
    for t in (0.0, 0.5):
        if spatial:
            add_spatial_step(domain, t, [
                grid(
                    topology("Triangle", 10 + r,
                             data_item([10 + r, 3], "out.h5:/e_p%d" % r,
                                       data_type="Int")),
                    geometry(data_item([8 + r, 3], "out.h5:/n_p%d" % r,
                                       precision=8)),
                    [attribute("u", data_item([8 + r], "out.h5:/u_p%d" % r,
                                              precision=8))],
                    name="subdomain_%d" % r,
                )
                for r in range(n_ranks)
            ])
        else:
            add_uniform_step(
                domain, t,
                topology("Tetrahedron", 12,
                         data_item([12, 4], "out.h5:/e", data_type="Int")),
                geometry(data_item([9, 3], "out.h5:/n", precision=8)),
                [attribute("u", data_item([9], "out.h5:/u", precision=8))],
            )
    return domain


@pytest.fixture
def archive(tmp_path):
    path = tmp_path / "run.ymf"
    write_ymf(make_domain(), path)
    return path


# --------------------------------------------------------------------------
# the conversion itself
# --------------------------------------------------------------------------


def test_the_output_defaults_to_the_archive_with_an_xmf_suffix(archive):
    destination = ymf2xdmf(archive)
    assert destination == archive.with_suffix(".xmf")
    assert destination.exists()


def test_an_explicit_output_path_is_honoured(archive, tmp_path):
    destination = ymf2xdmf(archive, tmp_path / "elsewhere.xmf")
    assert destination == tmp_path / "elsewhere.xmf"
    assert destination.exists()


def test_the_converted_document_is_xdmf(archive):
    root = parse(ymf2xdmf(archive)).getroot()
    assert root.tag == "Xdmf"
    collection = root.find("Domain").find("Grid")
    assert collection.attrib["CollectionType"] == "Temporal"
    assert len(collection.findall("Grid")) == 2


def test_the_conversion_preserves_every_dataset_reference(archive):
    """The .xmf must point at the same heavy data the .ymf does.

    Nothing is copied -- the .xmf references the same HDF5 datasets -- so a
    dropped or altered reference is the failure that yields a file which
    opens and shows the wrong thing.
    """
    domain, _ = read_ymf(archive)
    expected = []
    for collection in domain["TimeCollections"]:
        for step in collection["Data"]:
            for g in step.get("SpatialCollection", [step]):
                expected.append(g["Topology"]["DataItem"]["Data"])
                expected.append(g["Geometry"]["DataItem"]["Data"])
                for attr in g.get("Attributes", []):
                    expected.append(attr["DataItem"]["Data"])

    root = parse(ymf2xdmf(archive)).getroot()
    found = [(item.text or "").strip() for item in root.iter("DataItem")]
    assert found == expected


def test_a_per_subdomain_archive_converts_to_spatial_collections(tmp_path):
    path = tmp_path / "run.ymf"
    write_ymf(make_domain(spatial=True, n_ranks=3), path)
    root = parse(ymf2xdmf(path)).getroot()
    steps = root.find("Domain").find("Grid").findall("Grid")
    assert len(steps) == 2
    for step in steps:
        assert step.attrib["CollectionType"] == "Spatial"
        assert len(step.findall("Grid")) == 3


def test_the_extra_payload_survives_the_conversion(tmp_path):
    """Problem/vvuq metadata rides along in the <Information> element.

    So converting for a viewer does not discard the physics description,
    and read_xdmf can recover it.
    """
    from ymf.xdmf import read_xdmf

    extra = {"Problem": {"name": "poisson", "provenance": "proteus_derived"}}
    path = tmp_path / "run.ymf"
    write_ymf(make_domain(), path, extra=extra)
    _, recovered = read_xdmf(ymf2xdmf(path))
    assert recovered == extra


def test_an_empty_archive_is_refused(tmp_path):
    path = tmp_path / "empty.ymf"
    write_ymf({}, path, validate=False)
    with pytest.raises(YmfArchiveError, match="nothing to convert"):
        ymf2xdmf(path)


def test_a_malformed_archive_is_refused_by_default(tmp_path):
    domain = make_domain()
    domain["TimeCollections"][0]["Data"][0]["Attributes"][0]["Center"] = "Nowhere"
    path = tmp_path / "bad.ymf"
    write_ymf(domain, path, validate=False)
    with pytest.raises(YmfArchiveError, match="Center"):
        ymf2xdmf(path)


def test_validation_can_be_skipped_deliberately(tmp_path):
    domain = make_domain()
    domain["TimeCollections"][0]["Data"][0]["Attributes"][0]["Center"] = "Nowhere"
    path = tmp_path / "bad.ymf"
    write_ymf(domain, path, validate=False)
    assert ymf2xdmf(path, validate=False).exists()


# --------------------------------------------------------------------------
# the command line
# --------------------------------------------------------------------------


def test_the_cli_converts_and_reports_where(archive, capsys):
    assert ymf2xdmf_main([str(archive)]) == 0
    assert "wrote" in capsys.readouterr().out


def test_the_cli_honours_an_output_path(archive, tmp_path, capsys):
    destination = tmp_path / "out.xmf"
    assert ymf2xdmf_main([str(archive), "-o", str(destination)]) == 0
    assert destination.exists()


def test_verbose_describes_what_the_archive_holds(archive, capsys):
    assert ymf2xdmf_main([str(archive), "-v"]) == 0
    out = capsys.readouterr().out
    assert "Mesh Spatial_Domain" in out
    assert "2 steps" in out
    assert "Tetrahedron" in out
    assert "fields: u" in out


def test_a_missing_file_is_a_usage_error(tmp_path):
    with pytest.raises(SystemExit) as exc:
        ymf2xdmf_main([str(tmp_path / "nope.ymf")])
    assert exc.value.code == 2      # argparse's usage-error code


def test_a_malformed_archive_exits_nonzero_without_a_traceback(tmp_path, capsys):
    """A bad input file is the user's problem to see, not a stack trace."""
    domain = make_domain()
    domain["TimeCollections"][0]["Data"][0]["Attributes"][0]["Center"] = "Nowhere"
    path = tmp_path / "bad.ymf"
    write_ymf(domain, path, validate=False)
    assert ymf2xdmf_main([str(path)]) == 1
    assert "ymf2xdmf:" in capsys.readouterr().err


def test_the_cli_skips_validation_when_asked(tmp_path):
    domain = make_domain()
    domain["TimeCollections"][0]["Data"][0]["Attributes"][0]["Center"] = "Nowhere"
    path = tmp_path / "bad.ymf"
    write_ymf(domain, path, validate=False)
    assert ymf2xdmf_main([str(path), "--no-validate"]) == 0


# --------------------------------------------------------------------------
# the installed console script
# --------------------------------------------------------------------------


def test_the_installed_console_script_runs(archive):
    """Declared in pyproject as ymf2xdmf, so it must actually be callable."""
    import shutil

    executable = shutil.which("ymf2xdmf")
    if executable is None:
        pytest.skip("ymf2xdmf is not on PATH (package not installed)")
    result = subprocess.run([executable, str(archive)],
                            capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    assert archive.with_suffix(".xmf").exists()


def test_the_converter_needs_no_validation_extras(archive):
    """ymf2xdmf must work in a core-only install.

    It reads YAML and writes XML; requiring strictyaml or pint to convert
    would defeat the point of the core/front-end split.
    """
    code = (
        "import sys\n"
        "from importlib.abc import MetaPathFinder\n"
        "class Block(MetaPathFinder):\n"
        "    def find_spec(self, name, path=None, target=None):\n"
        "        if name.split('.')[0] in ('strictyaml', 'pint', 'sympy'):\n"
        "            raise ImportError(name)\n"
        "        return None\n"
        "sys.meta_path.insert(0, Block())\n"
        "from ymf.cli import ymf2xdmf\n"
        "print(ymf2xdmf(%r))\n" % (str(archive),)
    )
    result = subprocess.run([sys.executable, "-c", code],
                            capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
