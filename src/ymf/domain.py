"""The domain: the region a model is posed on.

A domain in YMF is more general than a mesh. XDMF's domain is discrete --
vertex coordinates (geometry) and connectivity (topology) -- and that is one
way to give a domain, not the only one. The same region can be stated as a
box, as a boundary representation, as an already-discretized mesh, or as a
reference to a file in some other language entirely: a ``.geo``, an STL, a
CAD model. Many domains are *technically* expressible in YAML and should
still be a file link, because inlining a CAD model helps nobody.

So a domain is a tagged union:

    kind: box                 lower/upper -- a parametric region
    kind: piecewise_linear    vertices, segments/facets -- a boundary
                              representation (Proteus's PSLG and PLC)
    kind: discrete            topology + geometry -- already a mesh, and
                              the one kind XDMF can hold directly
    kind: file                format + path -- someone else's language

``kind`` is written, not inferred. A reader should know what it is holding
before it reads it, and whoever writes a domain -- a person or a model --
should have to say what they are making.

The short form stays: ``domain: "Ω = [0, 1] × [0, 1]"`` is the notation a
spec is written in, and means ``kind: box``. ymf.symbolic parses it.

**What a domain is not.** It is not how to mesh that region, and it is not
the state of the object that happens to hold it. Proteus's Domain classes
carry both -- MeshOptions, he, use_gmsh alongside shape_list,
auxiliaryVariables, use_plex -- and neither belongs here: meshing policy
has a home already in ``discretizations[].mesh``, and runtime state is not
data. What is represented here is the region and how its boundary is
named, which is what the rest of a spec needs to refer to.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Set

__all__ = ["YmfDomainError", "DOMAIN_KINDS", "validate_domain_block", "box_bounds"]


class YmfDomainError(ValueError):
    """A domain block is not a well-formed domain of its stated kind."""


#: Keys every kind may carry.
_SHARED: Set[str] = {"kind", "dimension", "units", "boundary_tags",
                     "name", "time"}

#: What each kind requires, and what else it may hold.
DOMAIN_KINDS: Dict[str, Dict[str, Set[str]]] = {
    "box": {
        "required": {"lower", "upper"},
        "optional": set(),
    },
    "piecewise_linear": {
        "required": {"vertices"},
        # segments bound a 2-D region (PSLG); facets and volumes a 3-D one
        # (PLC). holes and regions mark what to leave out and what to tag.
        "optional": {"segments", "facets", "volumes", "holes", "regions",
                     "vertex_flags", "segment_flags", "facet_flags"},
    },
    "discrete": {
        "required": {"topology", "geometry"},
        "optional": set(),
    },
    "file": {
        "required": {"format", "path"},
        "optional": set(),
    },
}


def validate_domain_block(domain: Any, where: str = "domain") -> None:
    """Check a domain against its stated kind.

    A string is the notation short form and is left to ymf.symbolic. A
    mapping must say what kind it is, carry that kind's required keys, and
    carry nothing belonging to another kind -- a ``kind: box`` with
    ``vertices`` is a confusion worth catching rather than ignoring.
    """
    if isinstance(domain, str):
        return                      # the notation form; ymf.symbolic parses it
    if not isinstance(domain, dict):
        raise YmfDomainError(
            "%s: a domain is the notation string (\"Ω = [0, 1] × [0, 1]\") or a "
            "mapping with a kind; found %s" % (where, type(domain).__name__))

    kind = domain.get("kind")
    if kind is None:
        raise YmfDomainError(
            "%s: no kind. A domain says what it is: %s"
            % (where, ", ".join(sorted(DOMAIN_KINDS))))
    if kind not in DOMAIN_KINDS:
        raise YmfDomainError(
            "%s.kind: %r is not one of %s"
            % (where, kind, ", ".join(sorted(DOMAIN_KINDS))))

    spec = DOMAIN_KINDS[kind]
    missing = sorted(spec["required"] - set(domain))
    if missing:
        raise YmfDomainError(
            "%s: kind %r needs %s" % (where, kind, ", ".join(missing)))

    allowed = _SHARED | spec["required"] | spec["optional"]
    stray = sorted(set(domain) - allowed)
    if stray:
        elsewhere = {k: n for n, s in DOMAIN_KINDS.items()
                     for k in s["required"] | s["optional"] if k in stray}
        detail = "; ".join("%s belongs to kind %r" % (k, elsewhere[k])
                           for k in stray if k in elsewhere)
        raise YmfDomainError(
            "%s: kind %r does not take %s%s"
            % (where, kind, ", ".join(stray), ". " + detail if detail else ""))

    if kind == "box":
        lower, upper = domain["lower"], domain["upper"]
        if len(lower) != len(upper):
            raise YmfDomainError(
                "%s: lower has %d components and upper %d"
                % (where, len(lower), len(upper)))
        if not lower:
            raise YmfDomainError("%s: lower and upper are empty" % (where,))
        for i, (a, b) in enumerate(zip(lower, upper)):
            if not b > a:
                raise YmfDomainError(
                    "%s: upper[%d] = %r is not above lower[%d] = %r"
                    % (where, i, b, i, a))
        dim = domain.get("dimension")
        if dim is not None and dim != len(lower):
            raise YmfDomainError(
                "%s: dimension %r but %d components" % (where, dim, len(lower)))

    time = domain.get("time")
    if time is not None:
        if len(time) != 2 or not time[1] > time[0]:
            raise YmfDomainError(
                "%s.time: expected [t0, t1] with t1 above t0, got %r"
                % (where, time))

    if kind == "piecewise_linear" and not domain["vertices"]:
        raise YmfDomainError("%s: no vertices" % (where,))


def box_bounds(domain: Any) -> Optional[List[List[float]]]:
    """``[lower, upper]`` if this domain is a box in mapping form, else None.

    The notation short form returns None: it needs ymf.symbolic to parse,
    and this module stays dependency-free.
    """
    if isinstance(domain, dict) and domain.get("kind") == "box":
        return [list(domain["lower"]), list(domain["upper"])]
    return None
