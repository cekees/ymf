"""How a discretization reduces the model to something computable: its method.

A discretization entry names its ``method`` -- how the infinite-dimensional
problem becomes a finite, computable one -- and each method knows what an
entry needs before it can plan an output:

    mesh      finite elements (or volumes) on a mesh: ``mesh: {cells | h,
              levels}``. The only method implemented; the default, so
              every entry written before ``method`` existed is a mesh one.
    neural    a neural network trained on the equations (PINNs, weak
              adversarial networks): ``neural: {width, depth, points}``.
    particle  a meshless particle method (RKPM, SPH):
              ``particle: {spacing, kernel}``.

``neural`` and ``particle`` are schemas so far: their fields validate and
:func:`Method.missing` names what an entry lacks, but they plan no outputs
until their fields are settled together with a solver that runs them. A
spec using them can be written and reach ``well-posed``, not
``discretized``.

Each method also says which of its fields a study may list (each value one
realization, see ymf.closure), and how a realization reads in an output's
key. Whether a given *solver* can run a method is the solver's question,
not this module's: ``ymf_state`` must not depend on what is installed.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Dict, List, Tuple

__all__ = ["Missing", "Method", "METHODS", "DEFAULT", "method_of", "missing"]


@dataclass(frozen=True)
class Missing:
    """Something a discretization entry lacks before it can plan an output:
    ``what`` (a noun phrase, "a mesh") and ``hint`` (how to supply it)."""
    what: str
    hint: str


#: studied by every method that steps in time
_DT = ("dt", ("discretization", "dt"))


def _get(entry: Dict[str, Any], path: Tuple[str, ...]) -> Any:
    for part in path:
        if not isinstance(entry, dict) or part not in entry:
            return None
        entry = entry[part]
    return entry


@dataclass(frozen=True)
class Method:
    """One way of reducing the model to a computable problem."""
    name: str
    #: the entry's mapping that holds this method's parameters
    block: str
    #: (axis, path in the entry) that a study may list, in key order
    studied: Tuple[Tuple[str, Tuple[str, ...]], ...]
    #: values an axis takes when the block is present but the axis is not
    defaults: Dict[str, Any]
    #: what an entry of this method lacks (block-level checks are shared)
    lacks: Callable[[Dict[str, Any]], List[Missing]]
    #: False: the fields validate, but no output is planned yet
    implemented: bool = True

    def axes(self) -> Tuple[Tuple[str, Tuple[str, ...]], ...]:
        return self.studied + (_DT,)

    def missing(self, entry: Dict[str, Any]) -> List[Missing]:
        """What ``entry`` still lacks before it can plan an output."""
        out = []
        for other in METHODS.values():
            if other is not self and other.block in entry:
                out.append(Missing(
                    "a %s block that matches its method" % self.block,
                    "%s: belongs to method %s, and this entry's method is %s"
                    % (other.block, other.name, self.name)))
        out += self.lacks(entry)
        if not out and not self.implemented:
            out.append(Missing(
                "an implementation",
                "method %s is a schema only so far: its fields validate, but ymf "
                "plans no outputs for it yet" % (self.name,)))
        return out

    def key_parts(self, realization: Dict[str, Any]) -> List[str]:
        """How a realization reads in an output key, before the digest."""
        if self.name == "mesh":
            # the key format archives were first written with: cells or h,
            # then levels only when it is not 1
            size = "cells" if "cells" in realization else "h"
            parts = ["%s=%s" % (size, realization[size])]
            if realization.get("levels", 1) != 1:
                parts.append("levels=%s" % (realization["levels"],))
        else:
            parts = ["%s=%s" % (axis, realization[axis])
                     for axis, _ in self.studied if axis in realization]
        if "dt" in realization:
            parts.append("dt=%s" % (realization["dt"],))
        return parts


def _mesh_lacks(entry: Dict[str, Any]) -> List[Missing]:
    mesh = entry.get("mesh")
    if not isinstance(mesh, dict) or not ("cells" in mesh or "h" in mesh):
        return [Missing("a mesh", "give mesh: {h: ...} (or {cells: ...} on a box), "
                                  "or --h (--cells) on the command line")]
    if "cells" in mesh and "h" in mesh:
        return [Missing("one mesh size", "give mesh cells or h, not both: cells "
                                         "is for a box, h for any other region")]
    return []


def _fields_lacking(block: str, required: Tuple[str, ...], meaning: Dict[str, str]):
    def lacks(entry: Dict[str, Any]) -> List[Missing]:
        given = entry.get(block)
        given = given if isinstance(given, dict) else {}
        return [Missing("%s.%s" % (block, f), "give %s: {%s: ...}: %s" % (block, f, meaning[f]))
                for f in required if f not in given]
    return lacks


METHODS: Dict[str, Method] = {
    "mesh": Method(
        name="mesh", block="mesh",
        studied=(("cells", ("mesh", "cells")), ("h", ("mesh", "h")),
                 ("levels", ("mesh", "levels"))),
        defaults={"levels": 1},
        lacks=_mesh_lacks),
    "neural": Method(
        name="neural", block="neural",
        studied=(("width", ("neural", "width")), ("depth", ("neural", "depth")),
                 ("points", ("neural", "points"))),
        defaults={},
        lacks=_fields_lacking("neural", ("width", "depth", "points"), {
            "width": "neurons per hidden layer",
            "depth": "the number of hidden layers",
            "points": "the number of points the equations are enforced at"}),
        implemented=False),
    "particle": Method(
        name="particle", block="particle",
        studied=(("spacing", ("particle", "spacing")),),
        defaults={},
        lacks=_fields_lacking("particle", ("spacing", "kernel"), {
            "spacing": "the distance between particles",
            "kernel": "the kernel or shape function (e.g. cubic_spline, rk_linear)"}),
        implemented=False),
}

#: the method of an entry that names none: every entry written before
#: ``method`` existed reduces the model on a mesh
DEFAULT = "mesh"


def method_of(entry: Dict[str, Any]) -> Method:
    """The entry's method (``method:``, default mesh)."""
    name = entry.get("method", DEFAULT)
    if name not in METHODS:
        raise ValueError("discretization %r: method %r is not one of %s"
                         % (entry.get("name"), name, ", ".join(METHODS)))
    return METHODS[name]


def missing(entry: Dict[str, Any]) -> List[Missing]:
    """What ``entry`` lacks before it can plan an output (empty: nothing)."""
    return method_of(entry).missing(entry)
