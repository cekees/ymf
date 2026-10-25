"""Time outside, space inside: format 2 of the archive domain.

See docs/design-time-outside-space.md. The format 1 fixture was written by
ymf before the change (two spaces, three instants), with the .xmf it gave.
"""

from pathlib import Path

import pytest

from ymf.archive import (
    ARCHIVE_FORMAT_VERSION, YmfArchiveError, add_collection, add_instant,
    add_uniform_step, append_step, attribute, canonicalize_domain, data_item,
    geometry, new_domain, read_ymf, select_domain, series, space_names,
    spatial_grid, time_collections, topology, transpose_collections,
    uniform_grid, validate_domain, write_ymf,
)
from ymf.xdmf import build_xdmf_tree, parse_xdmf_domain, write_xdmf

FIXTURES = Path(__file__).parent / "fixtures"


def grid_at(name, t):
    return uniform_grid(
        topology("Triangle", 1, data_item([1, 3], "f.h5:/e")),
        geometry(data_item([3, 2], "f.h5:/n")),
        [attribute("u", data_item([3], "f.h5:/%s_%g" % (name, t)))])


def collection(name, times):
    return {"name": name, "data": [dict({"time": t}, **grid_at(name, t)) for t in times]}


def test_format_is_2():
    assert ARCHIVE_FORMAT_VERSION == 2


def test_equal_series_zip():
    times = transpose_collections([collection("A", [0.0, 0.5]), collection("B", [0.0, 0.5])])
    assert [i["time"] for i in times] == [0.0, 0.5]
    assert [list(i["spaces"]) for i in times] == [["A", "B"], ["A", "B"]]


def test_a_space_first_archived_late_is_absent_from_early_instants():
    times = transpose_collections([collection("A", [0.0, 0.5, 1.0]), collection("B", [1.0])])
    assert [list(i["spaces"]) for i in times] == [["A"], ["A"], ["A", "B"]]


def test_a_repeated_time_pairs_by_occurrence():
    times = transpose_collections([collection("A", [0.0, 1.0, 1.0]), collection("B", [0.0, 1.0, 1.0])])
    assert [i["time"] for i in times] == [0.0, 1.0, 1.0]
    assert all(list(i["spaces"]) == ["A", "B"] for i in times)


def test_close_but_unequal_times_are_distinct_instants():
    times = transpose_collections([collection("A", [0.1 + 0.2]), collection("B", [0.3])])
    assert len(times) == 2


def test_empty_and_single_step_collections():
    assert transpose_collections([collection("A", [])]) == []
    assert len(transpose_collections([collection("A", [2.0])])) == 1


def test_space_order_is_the_order_given():
    times = transpose_collections([collection("Z", [0.0]), collection("A", [0.0])])
    assert list(times[0]["spaces"]) == ["Z", "A"]


def test_the_format_1_fixture_reads_as_format_2():
    domain, _ = read_ymf(FIXTURES / "format1_two_spaces.ymf")
    assert list(domain) == ["times"]
    assert [i["time"] for i in domain["times"]] == [0.0, 0.5, 1.0]
    assert space_names(domain) == ["Mesh", "Mesh_c0p2"]


def test_format_1_converts_to_the_same_xmf_it_always_did(tmp_path):
    domain, _ = read_ymf(FIXTURES / "format1_two_spaces.ymf")
    write_xdmf(domain, tmp_path / "out.xmf")
    assert (tmp_path / "out.xmf").read_text() == (FIXTURES / "format1_two_spaces.xmf").read_text()


def test_xmf_round_trips_through_the_transpose(tmp_path):
    domain, _ = read_ymf(FIXTURES / "format1_two_spaces.ymf")
    write_xdmf(domain, tmp_path / "a.xmf")
    from ymf.xdmf import read_xdmf
    again, _ = read_xdmf(tmp_path / "a.xmf")
    assert again == canonicalize_domain(domain)


def test_write_then_read_is_format_2_and_equal(tmp_path):
    domain, _ = read_ymf(FIXTURES / "format1_two_spaces.ymf")
    write_ymf(domain, tmp_path / "x.ymf")
    assert (tmp_path / "x.ymf").read_text().startswith("ymf_archive_version: 2\n")
    again, _ = read_ymf(tmp_path / "x.ymf")
    assert again == canonicalize_domain(domain)


def test_the_format_1_api_still_builds_something_that_upgrades():
    d = new_domain("A")
    add_collection(d, "B")
    for t in (0.0, 1.0):
        for name in ("A", "B"):
            add_uniform_step(d, t, **{"topo": grid_at(name, t)["topology"], "geom": grid_at(name, t)["geometry"]},
                             collection=name)
    assert "time_collections" in d
    c = canonicalize_domain(d)
    assert [list(i["spaces"]) for i in c["times"]] == [["A", "B"], ["A", "B"]]


def test_format_2_steps_join_the_instant_at_their_time():
    d = new_domain()
    for t in (0.0, 1.0):
        for name in ("A", "B"):
            add_uniform_step(d, t, grid_at(name, t)["topology"], grid_at(name, t)["geometry"],
                             collection=name)
    assert [list(i["spaces"]) for i in d["times"]] == [["A", "B"], ["A", "B"]]


def test_an_instant_with_no_spaces_is_refused():
    with pytest.raises(YmfArchiveError, match=r"times\[0\]\.spaces"):
        canonicalize_domain({"times": [{"time": 0.0, "spaces": {}}]})


def test_errors_name_the_space_and_instant():
    bad = grid_at("A", 0.0)
    bad["topology"]["type"] = "Dodecahedron"
    d = new_domain()
    add_instant(d, 0.0, {"A": bad})
    with pytest.raises(YmfArchiveError, match=r"times\[0\]\.spaces\[A\]\.topology"):
        validate_domain(d)


def test_both_formats_in_one_domain_is_refused():
    with pytest.raises(YmfArchiveError, match="both"):
        canonicalize_domain({"times": [], "time_collections": []})


# --- append_step ---------------------------------------------------------


def header(tmp_path):
    path = tmp_path / "run.ymf"
    write_ymf(new_domain(), path)
    return path


def instant(t, *names):
    return {"time": t, "spaces": {n: grid_at(n, t) for n in names}}


def test_instants_append_and_read_back(tmp_path):
    path = header(tmp_path)
    for t in (0.0, 0.5):
        append_step(path, instant(t, "A", "B"))
    domain, _ = read_ymf(path)
    assert [i["time"] for i in domain["times"]] == [0.0, 0.5]
    assert space_names(domain) == ["A", "B"]


def test_a_space_may_first_appear_at_a_later_record(tmp_path):
    path = header(tmp_path)
    append_step(path, instant(0.0, "A"))
    append_step(path, instant(1.0, "A", "B"))
    domain, _ = read_ymf(path)
    assert [len(series(domain, n)) for n in ("A", "B")] == [2, 1]


def test_a_kill_loses_a_whole_instant_not_one_space_of_it(tmp_path):
    path = header(tmp_path)
    append_step(path, instant(0.0, "A", "B"))
    append_step(path, instant(1.0, "A", "B"))
    text = path.read_text()
    path.write_text(text[: text.rindex("B:")])           # cut inside the last instant
    with pytest.warns(RuntimeWarning):
        domain, _ = read_ymf(path)
    assert [i["time"] for i in domain["times"]] == [0.0]
    assert list(domain["times"][0]["spaces"]) == ["A", "B"]


def test_format_1_records_are_joined_into_instants(tmp_path):
    """What the Proteus writer does today: one record per collection per step."""
    path = header(tmp_path)
    for t in (0.0, 1.0):
        for name in ("A", "B"):
            append_step(path, name, dict({"time": t}, **grid_at(name, t)))
    domain, _ = read_ymf(path)
    assert [list(i["spaces"]) for i in domain["times"]] == [["A", "B"], ["A", "B"]]


def test_an_instant_record_in_a_format_1_document_is_refused(tmp_path):
    path = tmp_path / "old.ymf"
    path.write_text((FIXTURES / "format1_two_spaces.ymf").read_text())
    append_step(path, instant(2.0, "Mesh"))
    with pytest.raises(YmfArchiveError, match="format 1"):
        read_ymf(path)


def test_a_format_1_stream_in_flight_reads(tmp_path):
    path = tmp_path / "old.ymf"
    path.write_text("ymf_archive_version: 1\ndomain:\n  time_collections:\n  - name: A\n    data: []\n")
    for t in (0.0, 1.0):
        append_step(path, "A", dict({"time": t}, **grid_at("A", t)))
    domain, _ = read_ymf(path)
    assert [i["time"] for i in domain["times"]] == [0.0, 1.0]


def test_a_future_format_is_refused(tmp_path):
    path = tmp_path / "f.ymf"
    path.write_text("ymf_archive_version: 3\ndomain: {}\n")
    with pytest.raises(YmfArchiveError, match="version 3"):
        read_ymf(path)


# --- select and the views --------------------------------------------------


def three_by_two():
    d = new_domain()
    for t in (0.0, 1.0, 2.0):
        add_instant(d, t, {"A": grid_at("A", t), "B": grid_at("B", t)})
    return d


def test_steps_index_instants_across_spaces():
    picked = select_domain(three_by_two(), steps="-1")
    assert [i["time"] for i in picked["times"]] == [2.0]
    assert list(picked["times"][0]["spaces"]) == ["A", "B"]


def test_collections_filter_spaces_and_keep_archive_order():
    picked = select_domain(three_by_two(), collections=["B", "A"])
    assert space_names(picked) == ["A", "B"]
    only = select_domain(three_by_two(), collections=["B"])
    assert space_names(only) == ["B"] and len(only["times"]) == 3


def test_instants_without_a_chosen_space_are_dropped_before_steps_count():
    d = new_domain()
    add_instant(d, 0.0, {"A": grid_at("A", 0.0)})
    add_instant(d, 1.0, {"A": grid_at("A", 1.0), "B": grid_at("B", 1.0)})
    picked = select_domain(d, collections=["B"], steps="0")
    assert [i["time"] for i in picked["times"]] == [1.0]


def test_series_and_the_format_1_view_agree():
    d = three_by_two()
    assert [t for t, _ in series(d, "A")] == [0.0, 1.0, 2.0]
    view = time_collections(d)
    assert [c["name"] for c in view] == ["A", "B"]
    assert [s["time"] for s in view[1]["data"]] == [0.0, 1.0, 2.0]


def test_a_spatial_space_survives():
    d = new_domain()
    g = grid_at("A", 0.0)
    add_instant(d, 0.0, {"A": spatial_grid([dict(g, name="r0"), dict(g, name="r1")])})
    c = validate_domain(d)
    assert len(c["times"][0]["spaces"]["A"]["spatial_collection"]) == 2


def test_parse_xdmf_domain_returns_format_2():
    from xml.etree.ElementTree import ElementTree
    tree = build_xdmf_tree(three_by_two())
    domain = parse_xdmf_domain(tree.getroot())
    assert list(domain) == ["times"] and space_names(domain) == ["A", "B"]
