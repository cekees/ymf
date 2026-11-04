"""The reference graph: the edges between a document's branches.

A YMF document holds competing approaches to one problem and records how
they relate. Those relations are cross-references, and there are four:

    solution_paths.discretizations[].from_weak_form  -> weak_forms[].label
    vvuq.verification.reference.analytical           -> analytical[].name
    vvuq.verification.branches[]                     -> discretizations[].name
    vvuq.verification.expected_rates (keys)          -> discretizations[].name

Until now none of them was checked: a reference that named nothing still
validated, so a typo in a verification plan was simply ignored. They are
also where a document says what it *is* -- which weak form a
discretization implements is the difference between two branches of the
same problem and two different problems -- so leaving them unchecked made
the structure advisory.

Two of the four are single-valued, and those may be omitted when the graph
does not branch: a discretization in a document with one weak form can only
implement that one, and a verification plan in a document with one
analytical solution can only compare against that one. The omission is
filled in during composition, at the point the file is merged, and written
into the composed document, so the result has no implicit bindings and an
archive carries none either.

Filling in at that point rather than at the end is deliberate. A model with
one weak form is unambiguous when written; a problem extending it may add a
second. Resolving against the finished document would make the model's
omission retroactively ambiguous -- punishing a file for what a later file
did.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

__all__ = [
    "YmfReferenceError",
    "apply_defaults",
    "resolve",
    "weak_form_labels",
    "analytical_names",
    "discretization_names",
]


class YmfReferenceError(ValueError):
    """A reference names a branch the document does not have."""


def _problem(doc: Dict[str, Any]) -> Dict[str, Any]:
    return doc.get("model") or {}


def _paths(doc: Dict[str, Any]) -> Dict[str, Any]:
    return doc.get("solution_paths") or {}


def _verification(doc: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    block = (doc.get("vvuq") or {}).get("verification")
    return block if isinstance(block, dict) else None


def weak_form_labels(doc: Dict[str, Any]) -> List[str]:
    return [w["label"] for w in _problem(doc).get("weak_forms") or []
            if isinstance(w, dict) and "label" in w]


def analytical_names(doc: Dict[str, Any]) -> List[str]:
    return [a["name"] for a in _paths(doc).get("analytical") or []
            if isinstance(a, dict) and "name" in a]


def discretization_names(doc: Dict[str, Any]) -> List[str]:
    return [d["name"] for d in _paths(doc).get("discretizations") or []
            if isinstance(d, dict) and "name" in d]


def _only(candidates: List[str]) -> Optional[str]:
    """The single candidate, or None if there are none or several."""
    return candidates[0] if len(candidates) == 1 else None


# ---------------------------------------------------------------------------
# defaulting
# ---------------------------------------------------------------------------


def apply_defaults(doc: Dict[str, Any], source: str) -> List[Dict[str, Any]]:
    """Fill in omitted references where the graph does not branch.

    Returns one record per reference filled in, for ``composition.defaults``.
    Nothing is filled in where there is no candidate (a document may have
    discretizations and no weak forms) or more than one -- there, omitting
    the reference is a question only the author can answer, and
    :func:`resolve` asks it.
    """
    filled: List[Dict[str, Any]] = []

    label = _only(weak_form_labels(doc))
    if label is not None:
        for entry in _paths(doc).get("discretizations") or []:
            if isinstance(entry, dict) and "from_weak_form" not in entry:
                entry["from_weak_form"] = label
                filled.append({
                    "path": "solution_paths.discretizations[%s].from_weak_form"
                            % entry.get("name", "?"),
                    "value": label,
                    "because": "the document has one weak form",
                    "defaulted_by": source,
                })

    verification = _verification(doc)
    if verification is not None:
        reference = verification.get("reference")
        if isinstance(reference, dict) and "analytical" not in reference:
            name = _only(analytical_names(doc))
            if name is not None:
                reference["analytical"] = name
                filled.append({
                    "path": "vvuq.verification.reference.analytical",
                    "value": name,
                    "because": "the document has one analytical solution",
                    "defaulted_by": source,
                })
    return filled


# ---------------------------------------------------------------------------
# resolution
# ---------------------------------------------------------------------------


def _check(named: Any, candidates: List[str], path: str, what: str) -> None:
    if not isinstance(named, str):
        raise YmfReferenceError(
            "%s must name a %s; found %r" % (path, what, named))
    if named not in candidates:
        raise YmfReferenceError(
            "%s names %r, which is not %s %s in this document (have: %s)"
            % (path, named, "an" if what[:1] in "aeiou" else "a", what,
               ", ".join(candidates) or "none"))


def _check_equations(doc: Dict[str, Any]) -> None:
    """Named equations: each name once, and ``for`` naming the unknowns."""
    from ymf.normalize import equations_in_order, normalize_equations
    strong = _problem(doc).get("strong_form") or {}
    names = [e["name"] for e in normalize_equations(strong.get("equations") or [])
             if e["name"] is not None]
    twice = sorted({n for n in names if names.count(n) > 1})
    if twice:
        raise YmfReferenceError("model.strong_form.equations: %s named more than once"
                                % ", ".join(twice))
    if strong.get("equations") and strong.get("unknowns"):
        try:
            equations_in_order(strong, complete=doc.get("ymf_state") != "abstract")
        except ValueError as exc:
            raise YmfReferenceError("model.strong_form.%s" % exc) from exc


#: sharp functions of the notation a discretization may regularize
SHARP_FUNCTIONS = ("heaviside", "sign")


def _check_split(doc: Dict[str, Any], entry: Dict[str, Any]) -> None:
    """A discretization's own names: its sub-problems, bindings and
    regularizations, against the model and each other."""
    from ymf.normalize import normalize_equations, normalize_unknowns
    name = entry.get("name", "?")
    where = "solution_paths.discretizations[%s]" % (name,)
    subproblems = [s for s in entry.get("subproblems") or [] if isinstance(s, dict)]
    strong = _problem(doc).get("strong_form") or {}
    unknowns = [u["name"] for u in normalize_unknowns(strong.get("unknowns") or [])]
    coefficients = list(strong.get("coefficients") or {})
    model_names = set(unknowns) | set(coefficients)
    for symbol in entry.get("binds") or {}:
        if symbol in model_names:
            raise YmfReferenceError(
                "%s.binds.%s: %s is already a name in the model; a discretization "
                "binds symbols the model leaves to it" % (where, symbol, symbol))
    for function in entry.get("regularize") or {}:
        if function not in SHARP_FUNCTIONS:
            raise YmfReferenceError(
                "%s.regularize.%s: only the sharp functions %s are regularized"
                % (where, function, ", ".join(SHARP_FUNCTIONS)))
    if not subproblems:
        return
    names = [s["name"] for s in subproblems]
    twice = sorted({n for n in names if names.count(n) > 1})
    if twice:
        raise YmfReferenceError("%s.subproblems: %s named more than once"
                                % (where, ", ".join(twice)))
    equations = [e["name"] for e in normalize_equations(strong.get("equations") or [])]
    introduced: Dict[str, str] = {}
    for s in subproblems:
        new = (s.get("introduces") or {}).get("name")
        if new is None:
            continue
        if new in model_names or new in introduced:
            raise YmfReferenceError(
                "%s.subproblems[%s].introduces: %s is already a name (%s)"
                % (where, s["name"], new, "in the model" if new in model_names
                   else "introduced by %s" % introduced[new]))
        introduced[new] = s["name"]
    readable = model_names | set(introduced)
    solved_by: Dict[str, str] = {}
    for s in subproblems:
        here = "%s.subproblems[%s]" % (where, s["name"])
        if not (s.get("solves") or s.get("enforces") or s.get("introduces")
                or s.get("spin_up")):
            raise YmfReferenceError(
                "%s neither solves a model equation, enforces a constraint, "
                "introduces an unknown nor spins up: give one" % (here,))
        for equation in s.get("solves") or []:
            if None in equations:
                raise YmfReferenceError(
                    "%s.solves names %r, but the model's equations are not named; "
                    "name them ({name, for, equation}) to split them" % (here, equation))
            _check(equation, equations, here + ".solves", "model equation")
            if equation in solved_by:
                raise YmfReferenceError(
                    "equation %s is solved by both %s and %s in %s"
                    % (equation, solved_by[equation], s["name"], where))
            solved_by[equation] = s["name"]
        new = s.get("introduces") or {}
        for key, candidates in (("approximates", model_names), ("correction_of", unknowns)):
            if key in new:
                _check(new[key], sorted(candidates), "%s.introduces.%s" % (here, key),
                       "model unknown" if key == "correction_of" else
                       "model unknown or coefficient")
        for field, lag in (s.get("lags") or {}).items():
            copy_of = lag.get("of") if isinstance(lag, dict) else None
            if copy_of is not None:
                _check(copy_of, sorted(readable), "%s.lags.%s.of" % (here, field),
                       "model unknown, coefficient or introduced unknown")
                if field in readable:
                    raise YmfReferenceError(
                        "%s.lags.%s: a copy (of: %s) takes a name of its own, and %s is "
                        "already one" % (here, field, copy_of, field))
            else:
                _check(field, sorted(readable), "%s.lags" % (here,),
                       "model unknown, coefficient or introduced unknown")
            if isinstance(lag, dict) and "via" in lag:
                _check(lag["via"], names, "%s.lags.%s.via" % (here, field), "sub-problem")
        if "for" in s:
            _check(s["for"], sorted(set(unknowns) | set(introduced)), here + ".for",
                   "model or introduced unknown")
        for function in s.get("regularize") or {}:
            if function not in SHARP_FUNCTIONS:
                raise YmfReferenceError(
                    "%s.regularize.%s: only the sharp functions %s are regularized"
                    % (here, function, ", ".join(SHARP_FUNCTIONS)))
        for field in s.get("writes_back") or {}:
            _check(field, sorted(set(unknowns) | set(introduced)),
                   "%s.writes_back" % (here,), "model or introduced unknown")
        for symbol in s.get("binds") or {}:
            if symbol in readable:
                raise YmfReferenceError(
                    "%s.binds.%s: %s is already a name" % (here, symbol, symbol))
    if None not in equations:
        unsolved = [e for e in equations if e not in solved_by]
        if unsolved:
            raise YmfReferenceError(
                "%s: no sub-problem solves %s" % (where, ", ".join(unsolved)))


#: what each kind of condition needs: a value (value, formula or series),
#: none, or either
_NEEDS_VALUE = {"dirichlet": True, "neumann": True, "robin": True, "inflow": True,
                "characteristic": True, "periodic": False, "open": False}


def _check_conditions(doc: Dict[str, Any]) -> None:
    """Boundary conditions: their region, their variable, and what their kind needs."""
    strong = _problem(doc).get("strong_form") or {}
    unknowns = {u["name"]: int(u.get("rank") or 0)
                for u in _unknown_entries(strong.get("unknowns") or [])}
    coefficients = set(strong.get("coefficients") or {})
    regions = {r.get("name") for r in strong.get("boundary_regions") or []
               if isinstance(r, dict)}
    complete = doc.get("ymf_state") != "abstract"
    from ymf.normalize import boundary_conditions
    try:
        conditions = boundary_conditions(strong)
    except ValueError as exc:
        raise YmfReferenceError("model.strong_form.%s" % exc) from exc
    members = [dict(m, set=name) for name, entries in (strong.get("condition_sets") or {}).items()
               for m in entries or [] if isinstance(m, dict)]
    for bc in conditions + members:
        if not isinstance(bc, dict):
            continue
        where = "model.strong_form.%s[%s]" % (
            ("condition_sets.%s" % bc["set"]) if "region" not in bc else "boundary_conditions",
            ("variable=%r" % bc.get("variable")) if "region" not in bc else
            "region=%r, variable=%r%s" % (bc.get("region"), bc.get("variable"),
                                          ", set %s" % bc["set"] if bc.get("set") else ""))
        if "region" in bc and complete and regions and bc.get("region") not in regions:
            raise YmfReferenceError("%s: no boundary region named %r (have: %s)"
                                    % (where, bc.get("region"), ", ".join(sorted(regions))))
        variable, kind = bc.get("variable"), bc.get("type")
        if unknowns and variable not in unknowns and variable not in coefficients:
            raise YmfReferenceError("%s: %s is neither an unknown nor a coefficient"
                                    % (where, variable))
        given = any(bc.get(k) is not None for k in ("value", "formula", "series"))
        if kind == "slip" and variable in unknowns and unknowns[variable] != 1:
            raise YmfReferenceError("%s: slip gives the normal component of a vector, and "
                                    "%s is a scalar" % (where, variable))
        article = "an" if str(kind)[:1] in "aeiou" else "a"
        if _NEEDS_VALUE.get(kind) and not given:
            raise YmfReferenceError("%s: %s %s condition needs a value, formula or series"
                                    % (where, article, kind))
        if _NEEDS_VALUE.get(kind) is False and given:
            raise YmfReferenceError("%s: %s %s condition takes no value"
                                    % (where, article, kind))


def _unknown_entries(unknowns):
    from ymf.normalize import normalize_unknowns
    ranks = {u["name"]: u.get("rank") for u in unknowns if isinstance(u, dict)}
    return [dict(u, rank=ranks.get(u["name"])) for u in normalize_unknowns(unknowns)]


#: gauge kind -> (keys it needs, keys it may not have)
_GAUGE_KEYS = {
    "point": (("fields", "points"), ("lines", "tags", "quantities", "moment_point")),
    "line": (("fields", "lines"), ("points", "tags", "quantities", "moment_point")),
    "line_integral": (("fields", "lines"), ("points", "tags", "quantities", "moment_point")),
    "boundary_integral": (("tags", "quantities"), ("fields", "points", "lines")),
}


def _check_gauges(doc: Dict[str, Any]) -> None:
    """What each gauge's kind needs, and that its names exist."""
    model = _problem(doc)
    strong = model.get("strong_form") or {}
    names = ({u["name"] for u in _unknown_entries(strong.get("unknowns") or [])}
             | set(strong.get("coefficients") or {}))
    regions = {r.get("name") for r in strong.get("boundary_regions") or [] if isinstance(r, dict)}
    complete = doc.get("ymf_state") != "abstract"
    for i, gauge in enumerate(model.get("gauges") or []):
        if not isinstance(gauge, dict):
            continue
        where = "model.gauges[%s]" % (gauge.get("name") or i,)
        kind = gauge.get("kind")
        needs, refuses = _GAUGE_KEYS.get(kind, ((), ()))
        for key in needs:
            if not gauge.get(key):
                raise YmfReferenceError("%s: a %s gauge needs %s" % (where, kind, key))
        for key in refuses:
            if key in gauge:
                raise YmfReferenceError("%s: a %s gauge takes no %s" % (where, kind, key))
        if complete and names:
            for field in gauge.get("fields") or []:
                _check(field, sorted(names), where + ".fields", "unknown or coefficient")
        if complete and regions:
            for tag in gauge.get("tags") or []:
                _check(tag, sorted(regions), where + ".tags", "boundary region")
        points = list(gauge.get("points") or [])
        for line in gauge.get("lines") or []:
            if len(line) != 2:
                raise YmfReferenceError("%s: a line is its two ends, not %d points"
                                        % (where, len(line)))
            points.extend(line)
        if gauge.get("moment_point"):
            points.append(gauge["moment_point"])
        if len({len(p) for p in points}) > 1:
            raise YmfReferenceError("%s: its points do not all have the same dimension"
                                    % (where,))


def _check_regions(doc: Dict[str, Any]) -> None:
    """A boundary region is named by ``where`` or by ``tags``, exactly one; a
    tag is an int or a name the domain's ``boundary_tags`` resolves; a box's
    faces are named by position; a periodic region is found by position."""
    strong = _problem(doc).get("strong_form") or {}
    domain = strong.get("domain")
    is_box = isinstance(domain, str) or (isinstance(domain, dict)
                                         and domain.get("kind") == "box")
    names = (domain.get("boundary_tags") if isinstance(domain, dict) else None) or {}
    periodic = set()
    try:
        from ymf.normalize import boundary_conditions
        periodic = {bc.get("region") for bc in boundary_conditions(strong)
                    if bc.get("type") == "periodic"}
    except ValueError:
        pass                       # an undefined set: _check_conditions says so
    for region in strong.get("boundary_regions") or []:
        if not isinstance(region, dict):
            continue
        where = "model.strong_form.boundary_regions[%s]" % (region.get("name"),)
        given = [k for k in ("where", "tags") if region.get(k) is not None]
        if len(given) != 1:
            raise YmfReferenceError(
                "%s: give where (a predicate) or tags (the mesh's boundary tags), "
                "exactly one; this gives %s" % (where, " and ".join(given) or "neither"))
        if "tags" not in given:
            continue
        if is_box and domain is not None:
            raise YmfReferenceError(
                "%s: a box carries no boundary tags; name its faces with where "
                "(\"x = 0\")" % (where,))
        for tag in region["tags"]:
            if isinstance(tag, str) and not tag.lstrip("-").isdigit() and tag not in names:
                raise YmfReferenceError(
                    "%s: tag %r is not one of the domain's boundary_tags (have: %s)"
                    % (where, tag, ", ".join(sorted(names)) or "none"))
        if region.get("name") in periodic:
            raise YmfReferenceError(
                "%s: a periodic condition finds its faces by position, so the region "
                "needs where:, not tags:" % (where,))


def resolve(doc: Dict[str, Any]) -> None:
    """Check every reference, raising :class:`YmfReferenceError` on the first bad one."""
    _check_equations(doc)
    _check_regions(doc)
    _check_conditions(doc)
    _check_gauges(doc)
    for entry in _paths(doc).get("discretizations") or []:
        if isinstance(entry, dict):
            _check_split(doc, entry)
    labels = weak_form_labels(doc)
    discretizations = discretization_names(doc)

    for entry in _paths(doc).get("discretizations") or []:
        if not isinstance(entry, dict):
            continue
        named = entry.get("from_weak_form")
        where = "solution_paths.discretizations[%s].from_weak_form" % (
            entry.get("name", "?"),)
        if named is None:
            if labels:
                raise YmfReferenceError(
                    "%s is omitted, but this document has %d weak forms (%s), "
                    "so which one it implements has to be said"
                    % (where, len(labels), ", ".join(labels)))
            continue              # no weak forms at all: nothing to name
        _check(named, labels, where, "weak form label")

    verification = _verification(doc)
    if verification is None:
        return

    if "branch_ref" in (verification.get("reference") or {}):
        raise YmfReferenceError(
            "vvuq.verification.reference.branch_ref was renamed to "
            "vvuq.verification.reference.analytical")

    reference = verification.get("reference")
    if isinstance(reference, dict) and "analytical" in reference:
        _check(reference["analytical"], analytical_names(doc),
               "vvuq.verification.reference.analytical", "analytical solution")

    for i, named in enumerate(verification.get("branches") or []):
        _check(named, discretizations,
               "vvuq.verification.branches[%d]" % i, "discretization")

    for named in (verification.get("expected_rates") or {}):
        _check(named, discretizations,
               "vvuq.verification.expected_rates.%s" % named, "discretization")
