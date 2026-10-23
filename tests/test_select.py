"""Cherry-picking an archive: ymf2xmf --collection/--field/--steps, --show and -i."""
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
    select_domain,
    step_indices,
    topology,
    write_ymf,
)
from ymf.cli import describe, interactive_selection, ymf2xmf, ymf2xmf_main


def two_spaces(n_steps=4, spatial=False):
    """Two collections, several fields, n_steps instants; a Proteus-shaped archive."""
    domain = new_domain("Mesh Spatial_Domain")
    quadratic = {"name": "Mesh_c0p2_Lagrange", "data": []}
    domain["time_collections"].append(quadratic)
    for i in range(n_steps):
        t = 0.1 * i
        def fields(names):
            return [attribute(n, data_item([9], "out.h5:/%s_t%d" % (n, i), precision=8))
                    for n in names]
        if spatial:
            add_spatial_step(domain, t, [
                grid(topology("Triangle", 10, data_item([10, 3], "out.h5:/e_p%d" % r, data_type="Int")),
                     geometry(data_item([9, 3], "out.h5:/n_p%d" % r, precision=8)),
                     fields(["u", "p"]), name="p%d" % r)
                for r in range(2)], collection=0)
        else:
            add_uniform_step(domain, t,
                             topology("Triangle", 10, data_item([10, 3], "out.h5:/e", data_type="Int")),
                             geometry(data_item([9, 3], "out.h5:/n", precision=8)),
                             fields(["u", "p", "vof"]), collection=0)
        add_uniform_step(domain, t,
                         topology("Tri_6", 10, data_item([10, 6], "out.h5:/e2", data_type="Int")),
                         geometry(data_item([25, 3], "out.h5:/n2", precision=8)),
                         fields(["u2"]), collection=1)
    return domain


def names_of(domain):
    return [c["name"] for c in domain["time_collections"]]


def fields_of(step):
    return [a["name"] for a in step.get("attributes", [])]


# --------------------------------------------------------------------------
# step_indices
# --------------------------------------------------------------------------


@pytest.mark.parametrize("spec, expected", [
    ("-1", [3]), ("0", [0]), ("0,2", [0, 2]), ("1:3", [1, 2]), ("::2", [0, 2]),
    ("-2:", [2, 3]), ("2,0:2", [2, 0, 1]), ([3, 1], [3, 1]), (" 1 , 1 ", [1]),
])
def test_step_indices_read_like_python(spec, expected):
    assert step_indices(spec, 4) == expected


def test_a_step_the_archive_does_not_have_is_refused():
    with pytest.raises(YmfArchiveError, match="step 7 does not exist"):
        step_indices("7", 4)
    with pytest.raises(YmfArchiveError, match="neither an index nor a slice"):
        step_indices("last", 4)


# --------------------------------------------------------------------------
# select_domain
# --------------------------------------------------------------------------


def test_everything_is_kept_by_default():
    domain = two_spaces()
    assert select_domain(domain) == domain


def test_a_collection_is_picked_by_name_or_part_of_it():
    domain = two_spaces()
    assert names_of(select_domain(domain, collections=["Mesh Spatial_Domain"])) == ["Mesh Spatial_Domain"]
    assert names_of(select_domain(domain, collections=["c0p2"])) == ["Mesh_c0p2_Lagrange"]
    with pytest.raises(YmfArchiveError, match="no collection matches 'dgp1'"):
        select_domain(domain, collections=["dgp1"])


def test_fields_are_picked_on_every_grid():
    picked = select_domain(two_spaces(), fields=["p"])
    base, quadratic = picked["time_collections"]
    assert all(fields_of(step) == ["p"] for step in base["data"])
    # the quadratic space has no p: its grids keep their mesh and lose their fields
    assert all(fields_of(step) == [] for step in quadratic["data"])
    assert quadratic["data"][0]["topology"]["type"] == "Tri_6"


def test_fields_are_picked_inside_spatial_collections():
    picked = select_domain(two_spaces(spatial=True), fields=["u"], collections=["Spatial"])
    for step in picked["time_collections"][0]["data"]:
        assert [fields_of(g) for g in step["spatial_collection"]] == [["u"], ["u"]]


def test_an_unknown_field_is_refused_with_what_there_is():
    with pytest.raises(YmfArchiveError, match="no field named 'phi'.*'p', 'u', 'u2', 'vof'"):
        select_domain(two_spaces(), fields=["phi"])


def test_steps_apply_to_every_collection():
    picked = select_domain(two_spaces(), steps="-1")
    for c in picked["time_collections"]:
        assert [s["time"] for s in c["data"]] == [pytest.approx(0.3)]
    picked = select_domain(two_spaces(), steps="::2")
    assert [s["time"] for s in picked["time_collections"][0]["data"]] == [0.0, pytest.approx(0.2)]


def test_the_source_is_not_modified():
    domain = two_spaces()
    select_domain(domain, collections=["c0p2"], fields=["u2"], steps="0")
    assert names_of(domain) == ["Mesh Spatial_Domain", "Mesh_c0p2_Lagrange"]
    assert fields_of(domain["time_collections"][0]["data"][3]) == ["u", "p", "vof"]


# --------------------------------------------------------------------------
# ymf2xmf with a selection, --show, -i
# --------------------------------------------------------------------------


@pytest.fixture
def archive(tmp_path):
    path = tmp_path / "run.ymf"
    write_ymf(two_spaces(), path)
    return path


def grids_in(xmf):
    root = parse(xmf).getroot()
    return root.find("Domain").findall("Grid")


def test_the_conversion_can_cherry_pick(archive, tmp_path):
    out = ymf2xmf(archive, tmp_path / "part.xmf", collections=["Spatial"],
                  fields=["vof"], steps="-2:")
    (temporal,) = grids_in(out)
    assert temporal.attrib["Name"] == "Mesh Spatial_Domain"
    steps = temporal.findall("Grid")
    assert len(steps) == 2
    assert [a.attrib["Name"] for a in steps[0].findall("Attribute")] == ["vof"]
    # the archive itself is untouched
    domain, _ = read_ymf(archive)
    assert len(domain["time_collections"]) == 2


def test_the_flags_reach_the_conversion(archive, tmp_path, capsys):
    out = tmp_path / "part.xmf"
    assert ymf2xmf_main([str(archive), "-o", str(out), "-c", "c0p2", "-f", "u2", "-s", "0"]) == 0
    (temporal,) = grids_in(out)
    assert temporal.attrib["Name"] == "Mesh_c0p2_Lagrange"
    assert len(temporal.findall("Grid")) == 1


def test_a_bad_selection_is_an_error_not_a_traceback(archive, capsys):
    assert ymf2xmf_main([str(archive), "-f", "phi"]) == 1
    assert "no field named 'phi'" in capsys.readouterr().err


def test_show_describes_and_converts_nothing(archive, capsys):
    assert ymf2xmf_main([str(archive), "--show"]) == 0
    out = capsys.readouterr().out
    assert "1. Mesh Spatial_Domain: 4 steps, t=0..0.3, Triangle topology, fields: u, p, vof" in out
    assert "2. Mesh_c0p2_Lagrange" in out and "fields: u2" in out
    assert not archive.with_suffix(".xmf").exists()


def test_describe_lists_every_field_any_step_carries():
    domain = two_spaces(n_steps=2)
    # a field that appears only at the second step is still listed
    domain["time_collections"][0]["data"][1]["attributes"].append(
        attribute("late", data_item([9], "out.h5:/late", precision=8)))
    (base, _) = describe(domain)
    assert base.endswith("fields: u, p, vof, late")


def test_the_interactive_prompt_builds_the_selection(capsys):
    answers = iter(["1", "2 3", "-1"])
    collections, fields, steps = interactive_selection(two_spaces(), ask=lambda prompt: next(answers))
    assert (collections, fields, steps) == (["Mesh Spatial_Domain"], ["p", "vof"], "-1")
    out = capsys.readouterr().out
    assert "1. Mesh Spatial_Domain" in out and "1. u" in out


def test_blank_answers_keep_everything():
    answers = iter(["", "", ""])
    assert interactive_selection(two_spaces(), ask=lambda prompt: next(answers)) == (None, None, None)


def test_an_answer_that_is_not_a_choice_is_refused():
    with pytest.raises(YmfArchiveError, match="not one of 1..2"):
        interactive_selection(two_spaces(), ask=lambda prompt: "3")


def test_interactive_from_the_command_line(archive, tmp_path, monkeypatch):
    answers = iter(["2", "", ""])
    monkeypatch.setattr("builtins.input", lambda prompt="": next(answers))
    out = tmp_path / "part.xmf"
    assert ymf2xmf_main([str(archive), "-i", "-o", str(out)]) == 0
    (temporal,) = grids_in(out)
    assert temporal.attrib["Name"] == "Mesh_c0p2_Lagrange"
