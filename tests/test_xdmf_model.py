"""docs/xdmf-model.yaml describes the archive against XDMF; keep it true."""

from pathlib import Path

import pytest
import yaml

from ymf import archive

DOC = yaml.safe_load(
    (Path(__file__).parent.parent / "docs" / "xdmf-model.yaml").read_text(encoding="utf-8"))
ELEMENTS = DOC["XDMF_Elements"]
STATUSES = set(DOC["status_meanings"])


@pytest.mark.parametrize("where, values, constant", [
    (("Topology", "TopologyType"), None, archive.KNOWN_TOPOLOGY_TYPES),
    (("Geometry", "GeometryType"), None, archive.GEOMETRY_TYPES),
    (("Attribute", "AttributeType"), None, archive.ATTRIBUTE_TYPES),
    (("Attribute", "Center"), None, archive.CENTERINGS),
    (("XdmfItem", "NumberType"), None, archive.DATA_TYPES),
    (("XdmfItem", "Format"), None, archive.DATA_ITEM_FORMATS),
])
def test_closed_sets_match_the_archive_core(where, values, constant):
    element, attribute = where
    assert set(ELEMENTS[element]["attributes"][attribute]["values"]) == set(constant)


def test_sections_follow_the_xdmf_page_order():
    assert list(DOC["XML"]) == ["Elements", "XInclude", "XPath", "Entities"]
    assert list(ELEMENTS) == [
        "Xdmf", "Domain", "XdmfItem", "Grid", "Topology", "Geometry", "Attribute",
        "Set", "Time", "Information", "XML_Element_and_Default_XML_Attributes"]


def _statuses(node):
    if isinstance(node, dict):
        if "status" in node:
            yield node["status"]
        for value in node.values():
            yield from _statuses(value)


def test_every_status_is_one_of_the_defined_meanings():
    found = set(_statuses(DOC))
    assert found <= STATUSES, found - STATUSES
    assert "supported" in found and "not_supported" in found


def test_the_examples_are_valid_yaml():
    def examples(node):
        if isinstance(node, dict):
            for key, value in node.items():
                if key == "example":
                    yield value
                else:
                    yield from examples(value)
    found = list(examples(DOC))
    assert found
    for text in found:
        yaml.safe_load(text.replace("{...}", "{}").replace("[...]", "[]")
                       .replace("[step, step, ...]", "[]"))
