"""The two XDMF reads a Proteus archive reader needs, so Proteus imports no XML parser.

Proteus <= 1.9.x stored one ``<Grid>`` fragment per rank and step inside its
HDF5 files and wrote the collection order only into the ``.xmf`` beside them.
A current Proteus reads both when it hot starts from or converts such an
archive; these are the functions it reads them with.
"""
import pytest

from ymf.archive import YmfArchiveError
from ymf.xdmf import collection_names, read_grid_fragment

FRAGMENT = (
    '<Grid GridType="Uniform">'
    '<Time Value="2.5e-01" Name="3"/>'
    '<Topology Type="Triangle" NumberOfElements="2">'
    '<DataItem Format="HDF" DataType="Int" Precision="4" Dimensions="2 3">run.h5:/e3</DataItem>'
    '</Topology>'
    '<Geometry Type="XYZ">'
    '<DataItem Format="HDF" DataType="Float" Precision="8" Dimensions="4 3">run.h5:/n3</DataItem>'
    '</Geometry>'
    '<Attribute Name="u" AttributeType="Scalar" Center="Node">'
    '<DataItem Format="HDF" DataType="Float" Precision="8" Dimensions="4">run.h5:/u_t3</DataItem>'
    '</Attribute>'
    '</Grid>'
)


def test_a_fragment_with_a_time_gives_the_grid_and_the_time():
    grid, time = read_grid_fragment(FRAGMENT)
    assert time == 0.25
    assert grid["topology"]["data_item"]["data"] == "run.h5:/e3"
    assert grid["geometry"]["data_item"]["dimensions"] == [4, 3]
    assert [a["name"] for a in grid["attributes"]] == ["u"]
    assert "time" not in grid


def test_a_fragment_without_a_time_says_so():
    without = FRAGMENT.replace('<Time Value="2.5e-01" Name="3"/>', "")
    grid, time = read_grid_fragment(without)
    assert time is None
    assert grid["topology"]["type"] == "Triangle"


def test_bytes_are_read_as_well_as_text():
    # HDF5 hands the stored fragment back as bytes
    grid, time = read_grid_fragment(FRAGMENT.encode("utf-8"))
    assert time == 0.25 and grid["topology"]["number_of_elements"] == 2


def test_something_that_is_not_a_grid_is_refused():
    with pytest.raises(YmfArchiveError, match="expected a <Grid>"):
        read_grid_fragment("<Domain/>")
    with pytest.raises(YmfArchiveError, match="not XML"):
        read_grid_fragment("<Grid")


def test_collection_names_come_in_file_order(tmp_path):
    path = tmp_path / "run.xmf"
    path.write_text(
        '<?xml version="1.0" ?>\n'
        '<Xdmf Version="2.0"><Domain>'
        '<Grid Name="Mesh Spatial_Domain" GridType="Collection" CollectionType="Temporal"/>'
        '<Grid Name="EBMesh Spatial_Domain" GridType="Collection" CollectionType="Temporal"/>'
        '<Grid Name="Mesh_c0p2_Lagrange" GridType="Collection" CollectionType="Temporal">'
        '<Grid GridType="Collection" CollectionType="Spatial"><Grid GridType="Uniform"/></Grid>'
        '</Grid>'
        '</Domain></Xdmf>\n')
    # EBMesh sorts ahead of Mesh, which is why the order has to be read
    assert collection_names(path) == ["Mesh Spatial_Domain", "EBMesh Spatial_Domain",
                                      "Mesh_c0p2_Lagrange"]


def test_an_unnamed_collection_is_numbered(tmp_path):
    path = tmp_path / "run.xmf"
    path.write_text('<Xdmf><Domain><Grid GridType="Uniform"/>'
                    '<Grid Name="Mesh" GridType="Uniform"/></Domain></Xdmf>')
    assert collection_names(path) == ["Grid 0", "Mesh"]


def test_a_file_without_a_domain_has_no_collections(tmp_path):
    path = tmp_path / "run.xmf"
    path.write_text("<Xdmf/>")
    assert collection_names(path) == []
    path.write_text("<Xdmf")
    with pytest.raises(YmfArchiveError, match="not an XML document"):
        collection_names(path)
