"""A mesh's fingerprint: the digest a writer records in its provenance, so a reader can check it holds the
mesh it expects before reading data on it (C7.1 step I, decided by cekees 2026-11-05)."""

import pytest

from ymf.archive import (
    FINGERPRINT_METHOD,
    YmfArchiveError,
    provenance,
    read_ymf,
    validate_approximation,
    write_ymf,
)

np = pytest.importorskip("numpy")

from ymf.mesh import check_fingerprint, fingerprint, grid_fingerprint  # noqa: E402

from test_mesh import CELLS, NODES, mesh_domain, square_space  # noqa: E402

#: the digest of test_mesh's unit square; a change to the definition changes it, and must change the method
SQUARE_DIGEST = "4631cb7c6e4c040d31e317e86467b298eaea5d0c1cae3ea6e1e2bd8def99a91e"


def test_the_digest_is_of_the_numbers_in_their_order():
    digest = fingerprint(NODES, CELLS)
    assert len(digest) == 64 and digest == fingerprint(NODES[:, :2], CELLS.astype(np.int64))
    moved = NODES.copy()
    moved[2, 0] += 1e-15
    assert fingerprint(moved, CELLS) != digest
    assert fingerprint(NODES, CELLS[::-1]) != digest           # another cell order
    assert fingerprint(NODES, CELLS[:, [1, 2, 0]]) != digest   # another node order within a cell


def test_the_definition_is_pinned():
    assert fingerprint(NODES, CELLS) == SQUARE_DIGEST


def test_a_grid_records_its_fingerprint_and_reads_it_back(tmp_path):
    digest = fingerprint(NODES, CELLS)
    source = provenance("a mesh converter", fingerprint=digest.upper())
    assert source["fingerprint"] == {"method": FINGERPRINT_METHOD, "sha256": digest}
    domain = validate_approximation(mesh_domain(square_space(source=source)))
    write_ymf(domain, tmp_path / "square.ymf")
    back, _ = read_ymf(tmp_path / "square.ymf")
    space = validate_approximation(back)["times"][0]["spaces"]["Mesh Spatial_Domain"]
    assert grid_fingerprint(space) == digest
    assert check_fingerprint(space) == digest


def test_a_grid_holding_another_mesh_is_refused():
    space = square_space(source=provenance("t", fingerprint=fingerprint(NODES, CELLS[::-1])))
    with pytest.raises(YmfArchiveError, match="not the recorded"):
        check_fingerprint(space)
    # a writer that stored its rows in another order passes the digest of the input order
    assert check_fingerprint(space, digest=fingerprint(NODES, CELLS[::-1]))


def test_no_fingerprint_checks_nothing_and_another_method_is_refused():
    assert check_fingerprint(square_space(source=provenance("t"))) is None
    assert check_fingerprint(square_space()) is None
    space = square_space(source=provenance("t", fingerprint="ab" * 32))
    space["provenance"]["fingerprint"]["method"] = "some-other-1"
    with pytest.raises(YmfArchiveError, match="method 'some-other-1'"):
        check_fingerprint(space)


@pytest.mark.parametrize("block, match", [
    ({"tool": "t", "fingerprint": "ab" * 32}, r"a fingerprint is \{method, sha256\}"),
    ({"tool": "t", "fingerprint": {"method": FINGERPRINT_METHOD, "sha256": "12"}}, "not a SHA-256 digest"),
])
def test_bad_fingerprints_are_refused(block, match):
    with pytest.raises(YmfArchiveError, match=match):
        validate_approximation(mesh_domain(square_space(source=block)))
