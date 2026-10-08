"""A number in a spec is a number, whichever reader sees it.

Phase 1 of the semantics campaign. Three defects, one theme: a value
written in a spec used to arrive as a string, and which string depended on
which reader and how the number was spelled.

* ``coefficients`` was typed ``Any()``, under which strictyaml returns every
  scalar as a string. ``kappa: 1.0e-3`` loaded as ``'1.0e-3'``.
* PyYAML implements YAML 1.1, whose float grammar wants a decimal point and
  a signed exponent, so a spec read inside an archive turned ``1e-10`` into
  a string while the same spec read as a file did not.
* ``characteristic_scales`` was a closed set of six names, so a diffusivity
  or pressure scale had nowhere to go.
"""

import pytest

from ymf import load_ymf

strictyaml = pytest.importorskip("strictyaml")


SPEC = """\
kind: model
problem:
  name: "values"
  strong_form:
    provenance: human_edited
    equation_formulation: "−Δu = f"
    unknowns: ["u"]
    domain: "Ω"
    boundary_regions:
      - {name: "all", geometry: "∂Ω"}
    coefficients:
      plain: 0.5
      sci: 1.0e-3
      short_sci: 1e-10
      quoted: "1.0"
      formula: "π²"
      with_units: {value: 2.5, units: "m2/s"}
      declared_only: {units: "kg/m3"}
      as_formula: {formula: "(G, 0)"}
  characteristic_scales:
    length: {value: 1.0, units: "m"}
    diffusivity: {value: 1.0e-3, units: "m2/s"}
"""


@pytest.fixture
def spec(tmp_path):
    path = tmp_path / "values.ymf"
    path.write_text(SPEC, encoding="utf-8")
    return path


def coefficients(doc):
    return doc["problem"]["strong_form"]["coefficients"]


# --------------------------------------------------------------------------
# coefficients are typed
# --------------------------------------------------------------------------


@pytest.mark.parametrize("name,expected", [
    ("plain", 0.5),
    ("sci", 1.0e-3),
    ("short_sci", 1e-10),
    ("quoted", 1.0),
])
def test_a_numeric_coefficient_loads_as_a_number(spec, name, expected):
    value = coefficients(load_ymf(spec))[name]
    assert isinstance(value, float), "%s came back as %r" % (name, value)
    assert value == expected


def test_an_expression_coefficient_stays_a_string(spec):
    """Typing numbers must not take the formulas with it."""
    assert coefficients(load_ymf(spec))["formula"] == "π²"


def test_a_value_inside_a_map_is_a_number_too(spec):
    assert coefficients(load_ymf(spec))["with_units"] == {
        "value": 2.5, "units": "m2/s"}


def test_a_coefficient_may_declare_units_without_a_value(spec):
    """A model file names the slot; the problem extending it supplies one."""
    assert coefficients(load_ymf(spec))["declared_only"] == {"units": "kg/m3"}


def test_a_coefficient_may_be_a_formula_map(spec):
    assert coefficients(load_ymf(spec))["as_formula"] == {"formula": "(G, 0)"}


# --------------------------------------------------------------------------
# the two readers agree
# --------------------------------------------------------------------------


@pytest.mark.parametrize("text,expected", [
    ("1.0e-10", 1.0e-10),    # YAML 1.1 accepts this one
    ("1e-10", 1e-10),        # and these it does not
    ("1E-10", 1e-10),
    ("1e+10", 1e10),
    ("-1e-10", -1e-10),
    ("1.5e3", 1500.0),
])
def test_the_archive_reader_resolves_scientific_notation(text, expected):
    from ymf.archive import _load

    assert _load("v: %s" % text)["v"] == expected


def test_both_readers_agree_on_the_same_spec(spec):
    """The point of the exercise: one file, one meaning.

    strictyaml reads a spec written as a file; pyyaml reads the copy carried
    inside an archive. They used to disagree about what a number was.
    """
    from ymf.archive import _load

    by_pyyaml = _load(SPEC)["problem"]["strong_form"]["coefficients"]
    by_strictyaml = coefficients(load_ymf(spec))
    for name in ("plain", "sci", "short_sci"):
        assert by_pyyaml[name] == by_strictyaml[name]
        assert type(by_pyyaml[name]) is type(by_strictyaml[name]), name


def test_non_numbers_are_still_not_numbers():
    from ymf.archive import _load

    loaded = _load("a: notanumber\nb: 1.2.3\nc: 1e\nd: e10\n")
    assert all(isinstance(v, str) for v in loaded.values()), loaded


def test_pyyamls_own_loaders_are_left_alone():
    """The resolver goes on a subclass; other code in this process is not ours."""
    import yaml

    import ymf.archive  # noqa: F401 - registers the resolver

    assert yaml.safe_load("v: 1e-10")["v"] == "1e-10"


# --------------------------------------------------------------------------
# characteristic_scales is open
# --------------------------------------------------------------------------


def test_a_scale_outside_the_fixed_six_is_accepted(spec):
    scales = load_ymf(spec)["problem"]["characteristic_scales"]
    assert scales["diffusivity"] == {"value": 1.0e-3, "units": "m2/s"}


def test_an_extra_scale_resolves_under_its_own_name(tmp_path):
    """It is only useful if a formula can reach it."""
    pytest.importorskip("pint")
    from ymf.units import characteristic_scale_values

    path = tmp_path / "scales.ymf"
    path.write_text(SPEC, encoding="utf-8")
    values = characteristic_scale_values(load_ymf(path)["problem"])
    assert values["diffusivity"] == 1.0e-3
    assert values["length"] == 1.0
