"""The vocabulary of a finite element space: what a ``finite_element`` entry says beyond its family and order.

A space is stated as ``{family, order}`` plus, where they apply, the words below::

    fields: {family: CG, order: 3, basis: bernstein, bernstein_operator: false}
    fields: {family: CG, order: 4, basis: lagrange, nodes: gll}
    velocity: {family: RT, order: 0}
    geometry_order: 1                       # beside the spaces: the mesh map's polynomial degree

* ``family`` (schema): CG, DG, RT, BDM, RTF, N1curl.
* ``order``: the polynomial degree p of CG and DG (CG >= 1, DG >= 0). For the vector families it numbers the family's
  members, the lowest being RT0 (``order: 0``), BDM1 (``order: 1``) and N1curl (``order: 0``): ``RT, order 0`` is the
  lowest-order Raviart-Thomas space, whose normal trace on a facet is constant.
* ``basis`` (CG, DG): the polynomial space is P_p on simplices (Q_p on cubes); ``lagrange`` is the nodal basis,
  ``bernstein`` the Bernstein-Bezier basis (CG only: non-negative, sums to one, the dofs are control-point coefficients
  and not point values), ``monomial`` the DG basis of monomials. Default: ``lagrange`` for CG, ``monomial`` for DG.
* ``nodes`` (Lagrange basis): where the nodes sit, ``gll`` (warp-and-blend with Gauss-Lobatto-Legendre edge points on
  simplices, tensor GLL on cubes; well conditioned from p = 3) or ``equispaced``; they coincide for p <= 2. Default ``gll``.
* ``bernstein_operator`` (Bernstein basis): how data is put into the space. False (default) is the exact conversion of
  the degree-p interpolant (order p + 1); true takes control-point values (preserves bounds, second order).
* ``geometry_order`` (beside the spaces): the polynomial degree of the map from the reference cell to the element, 1 =
  affine (straight facets, constant Jacobian), 2 = quadratic, .. (iso- and super-parametric elements). Default 1.
* the vector families (RT, BDM, RTF, N1curl) have no ``basis`` or ``nodes``: their dofs are fixed by the family (normal
  flux densities for RT and BDM, tangential moments for N1curl), and their elements are mapped by a Piola transform
  (``contravariant_piola`` for RT, BDM, RTF; ``covariant_piola`` for N1curl; ``identity`` for CG and DG).

:func:`resolve` returns the space with every default filled in; two specs that resolve equal state the same space.
:data:`IMPLEMENTED_IN_PROTEUS` records what proteus's ``FemTools`` carries as of cekees/main 92f24ad9 (C8 rounds 1-5), as
information: ymf does not refuse a space proteus lacks, because a spec states a problem and a solver may not have caught up.
"""

from __future__ import annotations

from typing import Any, Dict

__all__ = ["SpaceError", "BASES", "NODES", "VECTOR_FAMILIES", "MAPS", "resolve", "check", "check_finite_element",
           "IMPLEMENTED_IN_PROTEUS", "implemented_in_proteus"]


class SpaceError(ValueError):
    """A space a spec states that does not make sense."""


BASES: Dict[str, str] = {
    "lagrange": "the nodal basis: dof i is the value at node i",
    "bernstein": "the Bernstein-Bezier basis: non-negative, sums to one; dofs are control-point coefficients",
    "monomial": "the basis of monomials (discontinuous spaces)",
}
NODES: Dict[str, str] = {
    "gll": "warp-and-blend nodes with Gauss-Lobatto-Legendre points on the edges (tensor GLL on cubes)",
    "equispaced": "equally spaced nodes (fine through p = 3; coincides with gll for p <= 2)",
}
#: family -> (lowest order, the Piola/pullback map of its elements, the bases it has)
_FAMILIES: Dict[str, Dict[str, Any]] = {
    "CG": {"min_order": 1, "map": "identity", "bases": ("lagrange", "bernstein"), "default_basis": "lagrange"},
    "DG": {"min_order": 0, "map": "identity", "bases": ("monomial", "lagrange"), "default_basis": "monomial"},
    "RT": {"min_order": 0, "map": "contravariant_piola", "bases": ()},
    "BDM": {"min_order": 1, "map": "contravariant_piola", "bases": ()},
    "RTF": {"min_order": 0, "map": "contravariant_piola", "bases": ()},
    "N1curl": {"min_order": 0, "map": "covariant_piola", "bases": ()},
}
VECTOR_FAMILIES = tuple(f for f, v in _FAMILIES.items() if not v["bases"])
MAPS = tuple(sorted({v["map"] for v in _FAMILIES.values()}))
_KEYS = {"family", "order", "basis", "nodes", "bernstein_operator"}

#: what proteus's FemTools carries (C8, cekees/main 92f24ad9): (family, basis) -> the orders. Information only.
IMPLEMENTED_IN_PROTEUS: Dict[Any, Any] = {
    ("CG", "lagrange"): "simplices: any p >= 1 (nodes gll or equispaced); cubes: p = 1, 2",
    ("CG", "bernstein"): "simplices: any p >= 1; cubes: p = 2",
    ("DG", "monomial"): "simplices: p = 0..6",
    ("DG", "lagrange"): "simplices: p = 1",
    ("RT", None): "simplices: order 0 (RT0)",
    ("BDM", None): "simplices: order 1 (BDM1)",
    "geometry_order": "1 (affine) only; higher waits for C7.1",
}


def resolve(space: Dict[str, Any], where: str = "finite_element") -> Dict[str, Any]:
    """The space with its defaults filled in: ``{family, order, basis, nodes, bernstein_operator, map}`` (``basis``,
    ``nodes``, ``bernstein_operator`` only where they apply). Raises :class:`SpaceError`."""
    unknown = sorted(set(space) - _KEYS)
    if unknown:
        raise SpaceError("%s: no key %s in a space; its keys are %s" % (where, ", ".join(unknown), ", ".join(sorted(_KEYS))))
    family = str(space["family"])
    spec = _FAMILIES[family]
    order = int(space["order"])
    if order < spec["min_order"]:
        raise SpaceError("%s: order %d: the lowest order of %s is %d" % (where, order, family, spec["min_order"]))
    out: Dict[str, Any] = {"family": family, "order": order, "map": spec["map"]}
    if not spec["bases"]:
        for key in ("basis", "nodes", "bernstein_operator"):
            if key in space:
                raise SpaceError("%s: %s has no %s: its dofs are fixed by the family (the %s map applies)"
                                 % (where, family, key, spec["map"]))
        return out
    basis = str(space.get("basis", spec["default_basis"]))
    if basis not in BASES:
        raise SpaceError("%s: basis %s is not one of %s" % (where, basis, ", ".join(BASES)))
    if basis not in spec["bases"]:
        raise SpaceError("%s: %s has no %s basis; its bases are %s" % (where, family, basis, ", ".join(spec["bases"])))
    out["basis"] = basis
    if "nodes" in space:
        nodes = str(space["nodes"])
        if nodes not in NODES:
            raise SpaceError("%s: nodes %s is not one of %s" % (where, nodes, ", ".join(NODES)))
        if basis != "lagrange":
            raise SpaceError("%s: nodes belong to a lagrange basis, and this one is %s" % (where, basis))
    if basis == "lagrange":
        out["nodes"] = str(space.get("nodes", "gll"))
    if "bernstein_operator" in space:
        if basis != "bernstein":
            raise SpaceError("%s: bernstein_operator belongs to a bernstein basis, and this one is %s" % (where, basis))
    if basis == "bernstein":
        out["bernstein_operator"] = bool(space.get("bernstein_operator", False))
    return out


def check(space: Dict[str, Any], where: str = "finite_element") -> Dict[str, Any]:
    """:func:`resolve`, under the name vocabulary checks use."""
    return resolve(space, where)


def check_finite_element(block: Any, where: str) -> None:
    """A ``finite_element`` block: each space it states, and its ``geometry_order``."""
    if not isinstance(block, dict):
        return
    for key, value in block.items():
        if key == "geometry_order":
            if int(value) < 1:
                raise SpaceError("%s.geometry_order: %s: the map's degree is at least 1 (1 = affine)" % (where, value))
        elif isinstance(value, dict):
            resolve(value, "%s.%s" % (where, key))


def implemented_in_proteus(space: Dict[str, Any]) -> bool:
    """Whether proteus's FemTools has this space on simplices (information, see :data:`IMPLEMENTED_IN_PROTEUS`)."""
    r = resolve(space)
    basis = r.get("basis")
    o = r["order"]
    f = r["family"]
    if f == "CG":
        return True
    if f == "DG":
        return (basis == "monomial" and o <= 6) or (basis == "lagrange" and o == 1)
    return (f == "RT" and o == 0) or (f == "BDM" and o == 1)
