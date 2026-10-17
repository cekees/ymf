"""The domain: the region a model is posed on.

Phase 5a of the semantics campaign. A domain used to be one thing -- the
notation string ``"Ω = [0, 1] × [0, 1]"`` -- which is a box and nothing
else. XDMF's domain is at the other extreme: discrete, vertex coordinates
and connectivity. Neither covers the middle, where most real domains live:
a boundary representation, or a reference to a ``.geo``, an STL or a CAD
model in a language YMF has no business reproducing.

So a domain is a tagged union, and ``kind`` is written rather than
inferred. What is represented is the region and how its boundary is named.
Not how to mesh it -- that is ``discretizations[].mesh`` -- and not the
state of whatever object happens to hold it.
"""

import pytest

from ymf import load_ymf
from ymf.domain import DOMAIN_KINDS, YmfDomainError, box_bounds, validate_domain_block

strictyaml = pytest.importorskip("strictyaml")


def spec(domain):
    return """\
ymf_state: abstract
model:
  name: "m"
  strong_form:
    provenance: llm_derived
    equation_formulation: "−Δu = f"
    unknowns: ["u"]
    domain: %s
    boundary_regions:
      - {name: "all", where: "∂Ω"}
solution_paths:
  analytical: []
  discretizations: []
""" % domain


def load(tmp_path, domain):
    path = tmp_path / "m.ymf"
    path.write_text(spec(domain), encoding="utf-8")
    return load_ymf(path)


# --------------------------------------------------------------------------
# the four kinds
# --------------------------------------------------------------------------


def test_the_notation_string_is_still_a_domain(tmp_path):
    """The short form stays: it is what a spec is written in."""
    doc = load(tmp_path, '"Ω = [0, 1] × [0, 1]"')
    assert doc["model"]["strong_form"]["domain"] == "Ω = [0, 1] × [0, 1]"


def test_a_box_in_mapping_form(tmp_path):
    doc = load(tmp_path, "{kind: box, lower: [0.0, 0.0], upper: [1.0, 2.0]}")
    assert box_bounds(doc["model"]["strong_form"]["domain"]) == [[0.0, 0.0], [1.0, 2.0]]


def test_a_boundary_representation(tmp_path):
    """Proteus's PSLG, which the notation form cannot say at all."""
    doc = load(tmp_path, """
      kind: piecewise_linear
      dimension: 2
      vertices: [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]
      segments: [[0, 1], [1, 2], [2, 3], [3, 0]]
      boundary_tags: {bottom: 1, right: 2, top: 3, left: 4}""")
    domain = doc["model"]["strong_form"]["domain"]
    assert len(domain["vertices"]) == 4
    assert domain["boundary_tags"]["bottom"] == 1


def test_a_file_reference_states_its_format(tmp_path):
    """Half of proteus.Domain is already a path; say so, and say in what."""
    doc = load(tmp_path, '{kind: file, format: "gmsh_geo", path: "flume.geo"}')
    domain = doc["model"]["strong_form"]["domain"]
    assert domain["format"] == "gmsh_geo" and domain["path"] == "flume.geo"


def test_a_discrete_domain_uses_the_archives_vocabulary(tmp_path):
    doc = load(tmp_path, """
      kind: discrete
      topology: {type: Triangle, number_of_elements: 2}
      geometry: {type: XYZ}""")
    assert doc["model"]["strong_form"]["domain"]["topology"]["type"] == "Triangle"


# --------------------------------------------------------------------------
# what it refuses
# --------------------------------------------------------------------------


def test_a_domain_must_say_what_kind_it_is():
    with pytest.raises(YmfDomainError, match="no kind"):
        validate_domain_block({"lower": [0], "upper": [1]})


def test_an_unknown_kind_lists_the_known_ones():
    with pytest.raises(YmfDomainError, match="piecewise_linear"):
        validate_domain_block({"kind": "mesh"})


def test_a_kind_missing_what_it_needs_is_told_what(tmp_path):
    with pytest.raises(YmfDomainError, match="kind 'file' needs format, path"):
        validate_domain_block({"kind": "file"})


def test_a_key_from_another_kind_is_refused_and_attributed():
    """A box with vertices is a confusion, not an extra."""
    with pytest.raises(YmfDomainError,
                       match=r"vertices belongs to kind 'piecewise_linear'"):
        validate_domain_block({"kind": "box", "lower": [0.0], "upper": [1.0],
                               "vertices": [[0.0]]})


@pytest.mark.parametrize("box,fragment", [
    ({"lower": [0.0, 0.0], "upper": [1.0]}, "2 components and upper 1"),
    ({"lower": [1.0], "upper": [0.0]}, "is not above"),
    ({"lower": [], "upper": []}, "empty"),
    ({"lower": [0.0], "upper": [1.0], "dimension": 3}, "dimension 3 but 1"),
])
def test_a_box_that_is_not_a_box(box, fragment):
    with pytest.raises(YmfDomainError, match=fragment):
        validate_domain_block(dict(box, kind="box"))


def test_a_time_interval_must_run_forwards():
    with pytest.raises(YmfDomainError, match=r"t1 above t0"):
        validate_domain_block({"kind": "box", "lower": [0.0], "upper": [1.0],
                               "time": [1.0, 0.0]})


def test_a_bad_domain_is_caught_when_the_document_loads(tmp_path):
    """Not only by a direct call: the check runs on every load path."""
    with pytest.raises(YmfDomainError, match="is not above"):
        load(tmp_path, "{kind: box, lower: [1.0], upper: [0.0]}")


# --------------------------------------------------------------------------
# running one
# --------------------------------------------------------------------------


def test_the_two_box_forms_give_the_same_problem():
    """The mapping form is a second spelling, not a second meaning."""
    pytest.importorskip("sympy")
    import copy

    from ymf.symbolic.problem import adr_problem

    doc = load_ymf("examples/poisson.ymf")
    mapped = copy.deepcopy(doc)
    mapped["model"]["strong_form"]["domain"] = {
        "kind": "box", "lower": [0.0, 0.0], "upper": [1.0, 1.0]}
    assert adr_problem(doc) == adr_problem(mapped)


def test_a_kind_nothing_solves_yet_says_so(tmp_path):
    """Validating and carrying a domain is not the same as running it."""
    pytest.importorskip("sympy")
    from ymf.symbolic.language import SymbolicError
    from ymf.symbolic.problem import adr_problem

    doc = load(tmp_path, '{kind: file, format: "stl", path: "hull.stl"}')
    doc["model"]["strong_form"]["equations"] = ["−Δu = 0"]
    with pytest.raises(SymbolicError, match="adr_problem can run a box"):
        adr_problem(doc)


def test_every_kind_is_covered_by_the_tests_above():
    """So a kind added later is not silently untested."""
    assert set(DOMAIN_KINDS) == {"box", "piecewise_linear", "discrete", "file"}
