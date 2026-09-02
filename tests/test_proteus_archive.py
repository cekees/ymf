"""Phase 0 acceptance: can YMF represent a real Proteus archive?

The fixtures in ``tests/fixtures/proteus/`` are **real, unmodified ``.xmf``
files produced by Proteus**, taken from a build tree and chosen to cover
every topology type Proteus was observed to emit (Triangle, Tri_6,
Polyline, Tetrahedron, Tet_10, Quadrilateral, Hexahedron).

The exit criterion for Phase 0 of the XML-to-YAML campaign is that a real
Proteus archive survives the whole pipeline the campaign introduces:

    Proteus .xmf  ->  YMF domain dict  ->  .ymf YAML  ->  domain dict
                  ->  .xmf  ->  domain dict

with the domain dict identical at every stage where it should be. If that
holds for real output, the model is wide enough to rebuild Proteus's writer
on top of. These tests are the thing that would fail if a later phase
narrowed the model.

One thing these fixtures caught that a hand-written fixture would not:
Proteus writes **one temporal collection per finite-element space**, not
one per archive. ``tri_6.xmf`` and ``tet_10.xmf`` each hold both
``Mesh Spatial_Domain`` (linear) and a ``Mesh_c0p2_Lagrange`` quadratic
space. 195 of the 814 archives surveyed had two collections; an earlier
version of the domain model had a single ``TimeCollection`` and silently
dropped the second. Hence ``TimeCollections`` being a list.

Note what is deliberately *not* asserted: byte-identity between Proteus's
``.xmf`` and YMF's. Proteus omits ``Precision`` on ``Int`` DataItems,
relying on XDMF's default of 4; YMF writes it explicitly. Those documents
are semantically identical and both valid, so the comparison is at the
level of the parsed model, not the bytes. See
:func:`test_reemitted_xmf_differs_from_proteus_only_in_explicit_defaults`,
which pins exactly that difference so it can't quietly become a bigger one.
"""

from pathlib import Path
from xml.etree.ElementTree import parse as et_parse

import pytest

from ymf.archive import (
    canonicalize_domain,
    read_ymf,
    validate_domain,
    write_ymf,
)
from ymf.xdmf import parse_xdmf_domain, write_xdmf

FIXTURE_DIR = Path(__file__).parent / "fixtures" / "proteus"
FIXTURES = sorted(FIXTURE_DIR.glob("*.xmf"))

# Every topology type observed across the 1121 .xmf files in the Proteus
# build tree these fixtures came from. Kept explicit so that if a fixture
# is ever dropped, the loss of coverage is visible rather than silent.
EXPECTED_TOPOLOGY_TYPES = {
    "Triangle",
    "Tri_6",
    "Polyline",
    "Tetrahedron",
    "Tet_10",
    "Quadrilateral",
    "Hexahedron",
}


def _parse(path):
    return parse_xdmf_domain(et_parse(path).getroot())


def test_fixtures_are_present():
    assert FIXTURES, f"no .xmf fixtures found in {FIXTURE_DIR}"


def test_fixtures_cover_every_observed_topology_type():
    found = set()
    for path in FIXTURES:
        for collection in _parse(path)["TimeCollections"]:
            for step in collection["Data"]:
                found.add(step["Topology"]["Type"])
    assert found == EXPECTED_TOPOLOGY_TYPES


@pytest.mark.parametrize("path", FIXTURES, ids=lambda p: p.stem)
def test_real_proteus_archive_parses_into_the_domain_model(path):
    domain = _parse(path)
    collections = domain["TimeCollections"]
    assert collections, "archive has no time collections"
    for tc in collections:
        assert tc["Name"], "time collection has no name"
        assert tc["Data"], "time collection has no timesteps"
        for step in tc["Data"]:
            assert "Time" in step
            assert step["Topology"]["NumberOfElements"] > 0
            assert step["Geometry"]["DataItem"]["Dimensions"]


@pytest.mark.parametrize("path", FIXTURES, ids=lambda p: p.stem)
def test_real_proteus_archive_passes_structural_validation(path):
    # If Proteus's own output fails YMF's structural checks, the checks are
    # wrong -- Proteus's archives are the ground truth here.
    validate_domain(_parse(path))


@pytest.mark.parametrize("path", FIXTURES, ids=lambda p: p.stem)
def test_real_proteus_archive_round_trips_through_yaml(path, tmp_path):
    domain = _parse(path)
    ymf_path = tmp_path / "archive.ymf"
    write_ymf(domain, ymf_path)
    assert read_ymf(ymf_path)[0] == canonicalize_domain(domain)


@pytest.mark.parametrize("path", FIXTURES, ids=lambda p: p.stem)
def test_real_proteus_archive_survives_the_full_pipeline(path, tmp_path):
    """xmf -> dict -> ymf -> dict -> xmf -> dict, unchanged throughout.

    This is the Phase 0 exit criterion.
    """
    original = canonicalize_domain(_parse(path))

    ymf_path = tmp_path / "archive.ymf"
    write_ymf(original, ymf_path)
    from_yaml, _ = read_ymf(ymf_path)
    assert from_yaml == original, "YAML round-trip changed the model"

    xmf_path = tmp_path / "reemitted.xmf"
    write_xdmf(from_yaml, xmf_path)
    from_xmf = _parse(xmf_path)
    assert from_xmf == original, "XDMF re-emission changed the model"


@pytest.mark.parametrize("path", FIXTURES, ids=lambda p: p.stem)
def test_reemitted_xmf_preserves_every_hdf5_dataset_reference(path, tmp_path):
    """No array reference may be lost, gained or altered.

    The failure mode this guards against is the one that produces an
    archive that opens without complaint and shows the wrong data: a
    DataItem whose reference got dropped or pointed at the wrong dataset.
    """

    def references(domain):
        refs = []
        for collection in domain["TimeCollections"]:
            refs.append(collection["Name"])
            for step in collection["Data"]:
                grids = step.get("SpatialCollection", [step])
                for g in grids:
                    refs.append(g["Topology"]["DataItem"].get("Data"))
                    refs.append(g["Geometry"]["DataItem"].get("Data"))
                    for attr in g.get("Attributes", []):
                        refs.append((attr["Name"], attr["DataItem"].get("Data")))
        return refs

    original = _parse(path)
    xmf_path = tmp_path / "reemitted.xmf"
    write_xdmf(canonicalize_domain(original), xmf_path)
    assert references(_parse(xmf_path)) == references(original)


@pytest.mark.parametrize("path", FIXTURES, ids=lambda p: p.stem)
def test_reemitted_xmf_differs_from_proteus_only_in_explicit_defaults(path, tmp_path):
    """Pin the one intentional difference from Proteus's own output.

    Proteus omits ``Precision`` on ``Int`` DataItems and lets XDMF's default
    of 4 apply; YMF always writes it. That is the *only* difference this
    test tolerates -- any other divergence in the DataItem attributes means
    the converter changed something it shouldn't have.
    """
    original_items = [
        elem.attrib for elem in et_parse(path).getroot().iter("DataItem")
    ]
    xmf_path = tmp_path / "reemitted.xmf"
    write_xdmf(canonicalize_domain(_parse(path)), xmf_path)
    reemitted_items = [
        elem.attrib for elem in et_parse(xmf_path).getroot().iter("DataItem")
    ]

    assert len(original_items) == len(reemitted_items)
    for original, reemitted in zip(original_items, reemitted_items):
        added = set(reemitted) - set(original)
        assert added <= {"Precision"}, f"unexpected added attributes: {added}"
        assert not set(original) - set(reemitted), "attributes were dropped"
        if "Precision" in added:
            assert original.get("DataType") == "Int"
            assert reemitted["Precision"] == "4"  # XDMF's default, made explicit
        for key in set(original) & set(reemitted):
            assert original[key] == reemitted[key], f"{key} changed"


def test_no_real_archive_uses_the_text_include_path():
    """Evidence, not just an assumption, about the --useTextArchive path.

    None of the 1121 ``.xmf`` files in the build tree these fixtures came
    from contained an ``xi:include``, and none used a Spatial collection.
    That is how the mis-targeted ``xi:include`` in Proteus's
    ``Archiver.py`` survived: nothing exercises that branch. Recorded here
    so the campaign's deprecation discussion rests on a measurement.
    """
    for path in FIXTURES:
        text = path.read_text()
        assert "xi:include" not in text
        assert 'CollectionType="Spatial"' not in text
