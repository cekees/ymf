"""Walk a YMF problem specification through the front-end.

    python examples/validate_spec.py                      # every example spec
    python examples/validate_spec.py examples/kovasznay_flow.ymf

For each document this prints what the validating front-end knows about
it: whether it passes the schema, its unknowns, the branch tree from
strong form to weak forms to discretizations, the unit check, and the
non-dimensionalization the units module derives. For a document built
with ``extends`` it also lists the files it was composed from and any
values a file overrode. Finally it breaks one document on purpose to show
what a unit error looks like.

Needs the front-end extras: ``pip install 'ymf[frontend]'``.
"""

from __future__ import annotations

import argparse
import copy
import sys
from pathlib import Path

import strictyaml

from ymf import check_units, load_ymf, non_dimensionalize, normalize_unknowns
from ymf.compose import YmfCompositionError

EXAMPLES = Path(__file__).resolve().parent


def describe(path: Path) -> bool:
    """Print a summary of one spec. Returns False if it fails validation."""
    print("=" * 72)
    print(path.relative_to(EXAMPLES) if path.is_relative_to(EXAMPLES) else path)
    print("=" * 72)
    try:
        doc = load_ymf(path)
    except (strictyaml.YAMLValidationError, YmfCompositionError) as exc:
        print("schema: INVALID\n%s" % (exc,))
        return False
    problem = doc["Problem"]
    if doc.get("kind") == "model":
        print("schema: valid model (not yet a problem: no domain, conditions or values)")
    else:
        print("schema: valid -- %s" % (problem["name"],))
    composition = doc.get("composition")
    if composition:
        print("composed from: %s" % " <- ".join(composition["sources"]))
        for o in composition.get("overrides", []):
            print("  override %s: %r -> %r (%s)" % (o["path"], o["was"], o["now"], o["set_by"]))

    strong = problem.get("strong_form", {})
    # An unknown without its own provenance inherits the strong form's.
    unknowns = normalize_unknowns(strong.get("unknowns", []), strong.get("provenance"))
    print("unknowns: %s" % ", ".join(
        "%s [%s] (%s)" % (u["name"], u.get("units") or "no units", u["provenance"])
        for u in unknowns))

    # The branch tree. Discretizations point at weak forms by label (or by
    # integer index -- both are accepted), and that link is what lets one
    # document hold several competing solution paths.
    print("branches:")
    print("  strong form: %s (%s)" % (strong.get("equation_formulation", "-"),
                                      strong.get("provenance", "-")))
    paths = doc.get("solution_paths", {"analytical": [], "discretizations": []})
    for index, wf in enumerate(problem.get("weak_forms", [])):
        print("    weak form %r (%s)" % (wf["label"], wf["provenance"]))
        for disc in paths["discretizations"]:
            # strictyaml hands an integer index back as a string
            if str(disc.get("from_weak_form")) in (wf["label"], str(index)):
                fe = disc.get("finite_element", {})
                spaces = ", ".join("%s %s%s" % (k, v["family"], v["order"])
                                   for k, v in fe.items())
                print("      discretization %r: %s; %s solver"
                      % (disc["name"], spaces or "-", disc["solver"].get("type", "-")))
    for exact in paths["analytical"]:
        print("    analytical %r: %s" % (exact["name"], exact["method"]))

    result = check_units(doc)
    print("unit check: %s" % ("ok" if result else "FAILED"))
    for issue in result:
        print("  %s" % (issue,))

    nd = non_dimensionalize(doc)
    if nd["characteristic_scales"]:
        print("characteristic scales: %s" % nd["characteristic_scales"])
    if nd["dimensionless_numbers"]:
        print("dimensionless numbers: %s" % nd["dimensionless_numbers"])
    if nd["substitutions"]:
        print("substitutions: %s" % nd["substitutions"])
    print()
    return True


def show_a_unit_error() -> None:
    """Break the heat-equation spec's boundary condition and re-check it."""
    print("=" * 72)
    print("a deliberate mistake: the wall temperature given in metres")
    print("=" * 72)
    doc = copy.deepcopy(load_ymf(EXAMPLES / "heat_equation.ymf"))
    doc["Problem"]["strong_form"]["boundary_conditions"][0]["units"] = "m"
    result = check_units(doc)
    print("unit check: %s" % ("ok" if result else "FAILED"))
    for issue in result:
        print("  %s" % (issue,))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("specs", nargs="*", type=Path,
                        help="YMF documents to check (default: every example)")
    args = parser.parse_args(argv)
    specs = [p.resolve() for p in args.specs] or sorted(EXAMPLES.rglob("*.ymf"))
    ok = all([describe(p) for p in specs])
    if not args.specs:
        show_a_unit_error()
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
