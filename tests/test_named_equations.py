"""Named equations: ``{name, for, equation}`` beside the bare strings.

A name is what a discretization's sub-problem refers to and what an
extending file merges onto; ``for`` says which unknown an equation
determines, so the order they are written in stops mattering.
"""

from pathlib import Path

import pytest

pytest.importorskip("strictyaml")

from ymf.normalize import equations_in_order, normalize_equations  # noqa: E402
from ymf.references import YmfReferenceError  # noqa: E402
from ymf.schema import load_ymf  # noqa: E402

EXAMPLES = Path(__file__).parent.parent / "examples"
NS = EXAMPLES / "navier_stokes"

MOMENTUM = "ρ ∂v/∂t + ∇·(ρ v⊗v) − μΔv + ∇p = f  in Ω"
CONTINUITY = "∇·v = 0  in Ω"


def _named(**extra):
    return [{"name": "continuity", "for": "p", "equation": CONTINUITY},
            {"name": "momentum", "for": "v", "equation": MOMENTUM, **extra}]


# --- the helpers, no schema -------------------------------------------------


def test_bare_strings_keep_their_order():
    strong = {"unknowns": ["v", "p"], "equations": [MOMENTUM, CONTINUITY]}
    assert equations_in_order(strong) == [MOMENTUM, CONTINUITY]


def test_for_puts_them_in_the_unknowns_order():
    strong = {"unknowns": ["v", "p"], "equations": _named()}
    assert equations_in_order(strong) == [MOMENTUM, CONTINUITY]


def test_a_named_equation_without_for_keeps_its_place():
    strong = {"unknowns": ["u"], "equations": [{"name": "heat", "equation": "-Δu = f"}]}
    assert normalize_equations(strong["equations"]) == [
        {"name": "heat", "for": None, "equation": "-Δu = f"}]
    assert equations_in_order(strong) == ["-Δu = f"]


@pytest.mark.parametrize("equations, message", [
    ([{"name": "momentum", "for": "v", "equation": MOMENTUM}, CONTINUITY],
     "give 'for' on every equation or on none"),
    ([{"name": "momentum", "for": "w", "equation": MOMENTUM},
      {"name": "continuity", "for": "p", "equation": CONTINUITY}],
     "is for 'w', which is not an unknown"),
    ([{"name": "a", "for": "v", "equation": MOMENTUM},
      {"name": "b", "for": "v", "equation": CONTINUITY}],
     "a and b are both for v"),
    ([{"name": "momentum", "for": "v", "equation": MOMENTUM}],
     "none is for p"),
])
def test_what_does_not_say_which_equation_is_for_which_unknown(equations, message):
    with pytest.raises(ValueError, match=message):
        equations_in_order({"unknowns": ["v", "p"], "equations": equations})


def test_an_incomplete_set_is_fine_where_the_chain_completes_it():
    strong = {"unknowns": ["v", "p"],
              "equations": [{"name": "momentum", "for": "v", "equation": MOMENTUM}]}
    assert equations_in_order(strong, complete=False) == [MOMENTUM]


# --- through the loader -----------------------------------------------------


def _poiseuille_with(tmp_path, equations_yaml, model_extra=""):
    """plane_poiseuille, extending a copy of the model whose equations are given."""
    model = (NS / "navier_stokes_model.ymf").read_text(encoding="utf-8")
    start = model.index("    equations:")
    end = model.index("    # The coefficients the equations use")
    (tmp_path / "navier_stokes_model.ymf").write_text(
        model[:start] + equations_yaml + model[end:] + model_extra, encoding="utf-8")
    (tmp_path / "plane_poiseuille.ymf").write_text(
        (NS / "plane_poiseuille.ymf").read_text(encoding="utf-8"), encoding="utf-8")
    return tmp_path / "plane_poiseuille.ymf"


NAMED_YAML = (
    "    equations:\n"
    "      - name: continuity\n"
    "        for: p\n"
    '        equation: "∇·v = 0  in Ω"\n'
    "      - name: momentum\n"
    "        for: v\n"
    '        equation: "ρ ∂v/∂t + ∇·(ρ v⊗v) − μΔv + ∇p = f  in Ω"\n')


def test_named_and_reordered_equations_give_the_same_problem(tmp_path):
    pytest.importorskip("sympy")
    from ymf.symbolic import adr_problem
    named = load_ymf(_poiseuille_with(tmp_path, NAMED_YAML))
    plain = load_ymf(NS / "plane_poiseuille.ymf")
    assert adr_problem(named)["adr"] == adr_problem(plain)["adr"]


def test_an_extending_file_replaces_an_equation_by_its_name(tmp_path):
    path = _poiseuille_with(tmp_path, NAMED_YAML)
    (tmp_path / "variant.ymf").write_text(
        "extends: plane_poiseuille.ymf\n"
        "model:\n"
        "  strong_form:\n"
        "    equations:\n"
        "      - name: momentum\n"
        '        equation: "ρ ∂v/∂t + ∇·(ρ v⊗v) − ∇·(μ∇v) + ∇p = f  in Ω"\n',
        encoding="utf-8")
    assert path.exists()
    doc = load_ymf(tmp_path / "variant.ymf")
    equations = doc["model"]["strong_form"]["equations"]
    assert [e["name"] for e in equations] == ["continuity", "momentum"]
    assert equations[1]["for"] == "v"                       # kept from the model
    assert "∇·(μ∇v)" in equations[1]["equation"]
    assert any(o["path"] == "model.strong_form.equations[momentum].equation"
               for o in doc["composition"]["overrides"])


def test_an_abstract_model_may_leave_an_unknown_to_the_files_that_extend_it(tmp_path):
    one = ("    equations:\n"
           "      - name: momentum\n"
           "        for: v\n"
           '        equation: "ρ ∂v/∂t + ∇·(ρ v⊗v) − μΔv + ∇p = f  in Ω"\n')
    path = _poiseuille_with(tmp_path, one)
    load_ymf(tmp_path / "navier_stokes_model.ymf")          # abstract: fine
    with pytest.raises(YmfReferenceError, match="none is for p"):
        load_ymf(path)                                     # well-posed: not


@pytest.mark.parametrize("equations_yaml, message", [
    (NAMED_YAML.replace("for: p", "for: q"), "is for 'q', which is not an unknown"),
    (NAMED_YAML.replace("name: continuity", "name: momentum"), "momentum named more than once"),
])
def test_bad_names_are_refused_on_load(tmp_path, equations_yaml, message):
    with pytest.raises(YmfReferenceError, match=message):
        load_ymf(_poiseuille_with(tmp_path, equations_yaml))


def test_the_flux_units_check_reads_named_equations(tmp_path):
    pytest.importorskip("sympy")
    pytest.importorskip("pint")
    from ymf import check_units
    named = load_ymf(_poiseuille_with(tmp_path, NAMED_YAML))
    plain = load_ymf(NS / "plane_poiseuille.ymf")
    assert ([str(i) for i in check_units(named, missing_units="warn").issues]
            == [str(i) for i in check_units(plain, missing_units="warn").issues])
