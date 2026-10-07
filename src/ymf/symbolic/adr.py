"""Strong form -> advection-diffusion-reaction (ADR) form.

Proteus, like many transport codes, solves systems in the form

    dm_i/dt + div(f_i - sum_k a_ik grad(phi_k)) + r_i + H_i = 0

one scalar equation ``i`` per component, with the mass ``m``, advective
flux ``f``, diffusion tensors ``a`` acting on potentials ``phi``, reaction
``r`` and Hamiltonian ``H`` (the part that is not in divergence form) all
functions of the unknowns, the coordinates and time. This module reads a
scalarized strong form (see :mod:`ymf.symbolic.language`) and sorts every
term into one of those slots:

=================================  ===================================
term                               slot
=================================  ===================================
``c * Dt(e)``                      mass, ``m += c e``
``c * D(e, j)``, e first order     advective flux, ``f_j += c e``
``c * D(b * D(phi, k), j)``        diffusion, ``a_jk -= c b``, potential phi
no derivatives                     reaction, ``r += term``
anything else first order          Hamiltonian, ``H += term``
=================================  ===================================

where ``c`` is a number. That table is ibvp's Proteus classifier
(https://github.com/ibvp/ibvp, MIT, Kloeckner and Kirby) on sympy.

It also works out what Proteus needs to assemble a Newton Jacobian: the
derivative of every coefficient with respect to every unknown, and whether
each dependence is constant, linear or nonlinear.

The result, :func:`adr_form`, is plain data -- numpy code strings and
flags, no sympy objects -- so a solver can consume it without importing
sympy. That is the same trust boundary the archive core keeps: the
front-end does the symbolic work once; the solver evaluates strings.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence, Tuple

import sympy
from sympy.printing.numpy import NumPyPrinter

from ymf.symbolic.language import D, Dt, Space, SymbolicError

__all__ = ["classify", "adr_form", "to_code", "grad_name"]


def grad_name(component: str, axis: int) -> str:
    """Name of the gradient component d(component)/dx_axis in emitted code."""
    return "grad_%s_%d" % (component, axis)


class _Printer(NumPyPrinter):
    """numpy code, with every Symbol printed as its plain name."""

    def __init__(self):
        super().__init__({"fully_qualified_modules": True})

    def _print_Symbol(self, expr):
        return expr.name

    def _print_Float(self, expr):
        # repr is the shortest text that reads back to the same double
        return repr(float(expr))


def to_code(expr) -> str:
    return _Printer().doprint(sympy.sympify(expr))


class Equation:
    """One classified scalar equation, still in sympy."""

    def __init__(self, dim: int):
        self.mass = sympy.Integer(0)
        self.advection = [sympy.Integer(0)] * dim
        #: potential component -> dim x dim tensor
        self.diffusion: Dict[str, sympy.Matrix] = {}
        self.reaction = sympy.Integer(0)
        self.hamiltonian = sympy.Integer(0)


def _split(term) -> Tuple[sympy.Expr, List[sympy.Expr]]:
    """(numeric factor, other factors) of a product."""
    factors = sympy.Mul.make_args(term)
    number = sympy.Mul(*[f for f in factors if f.is_number])
    return number, [f for f in factors if not f.is_number]


def _derivative_order(expr) -> int:
    if isinstance(expr, D):
        return 1 + _derivative_order(expr.args[0])
    if isinstance(expr, Dt):
        return _derivative_order(expr.args[0])
    return max((_derivative_order(a) for a in expr.args), default=0)


def _diffusion_terms(arg, space: Space, where: str):
    """Read ``sum_k b_k D(phi_k, k)`` -> [(b_k, phi_k, k)]."""
    out = []
    for term in sympy.Add.make_args(sympy.expand(arg, deep=False)):
        inner = [f for f in sympy.Mul.make_args(term) if isinstance(f, D)]
        rest = sympy.Mul(*[f for f in sympy.Mul.make_args(term) if not isinstance(f, D)])
        if len(inner) != 1 or rest.has(D, Dt):
            raise SymbolicError(
                "%s: %s is not a diffusive flux b*grad(phi) -- second derivatives "
                "must appear as div(b*grad(phi))" % (where, term))
        phi, k = inner[0].args
        if phi not in space.component_symbols:
            raise SymbolicError(
                "%s: the potential %s is not an unknown; only phi = an unknown is "
                "supported so far" % (where, phi))
        out.append((rest, phi, int(k)))
    return out


def classify(residual, space: Space, where: str = "equation",
             hamiltonian_gradients=()) -> Equation:
    """Sort the terms of one scalar residual into ADR slots.

    ``hamiltonian_gradients`` names unknowns whose bare gradient terms
    (``c * D(p, j)``) go to the Hamiltonian instead of the advective flux:
    ``∇p`` as ``H = ∂p/∂x_j``, not integrated by parts, rather than as the
    flux ``div(p I)``. The two are the same strong form and different weak
    forms.
    """
    eq = Equation(space.dim)
    allowed = set(space.component_symbols) | set(space.x) | {space.t}
    stray = sympy.sympify(residual).free_symbols - allowed
    if stray:
        raise SymbolicError("%s: undefined names %s (not a coordinate, an unknown "
                            "or a coefficient)" % (where, ", ".join(sorted(map(str, stray)))))
    for term in sympy.Add.make_args(sympy.expand(residual, deep=False)):
        number, factors = _split(term)
        ops = [f for f in factors if f.has(D, Dt)]
        rest = sympy.Mul(*[f for f in factors if not f.has(D, Dt)])
        if not ops:
            eq.reaction += term
            continue
        if len(ops) == 1 and rest == 1:
            op = ops[0]
            if isinstance(op, Dt):
                eq.mass += number * op.args[0]
                continue
            if isinstance(op, D):
                arg, axis = op.args[0], int(op.args[1])
                if isinstance(arg, sympy.Symbol) and arg.name in hamiltonian_gradients:
                    eq.hamiltonian += term
                    continue
                if not arg.has(D, Dt):
                    eq.advection[axis] += number * arg
                    continue
                for b, phi, k in _diffusion_terms(arg, space, where):
                    tensor = eq.diffusion.setdefault(
                        phi.name, sympy.zeros(space.dim, space.dim))
                    tensor[axis, k] += -number * b
                continue
        # Not in divergence form: must be first order, no time derivative.
        if term.has(Dt):
            raise SymbolicError("%s: time derivative in %s is not of the form c*dt(e)"
                                % (where, term))
        if _derivative_order(term) > 1:
            raise SymbolicError("%s: second derivatives in %s are not in "
                                "divergence form" % (where, term))
        for d in term.atoms(D):
            if d.args[0] not in space.component_symbols:
                raise SymbolicError("%s: in %s, only derivatives of an unknown may "
                                    "appear outside a divergence" % (where, term))
        eq.hamiltonian += term
    return eq


# --------------------------------------------------------------------------
# dependence and derivatives
# --------------------------------------------------------------------------


def _gradient_symbols(space: Space):
    """Replace D(u, k) by a plain symbol, for differentiating H."""
    subs = {}
    for comp in space.components:
        for k in range(space.dim):
            subs[D(sympy.Symbol(comp, real=True), k)] = sympy.Symbol(grad_name(comp, k), real=True)
    return subs


def _flag(expr, space: Space) -> Dict[str, str]:
    """{component: 'linear' | 'nonlinear'} for each unknown expr depends on.

    Linear means the derivative is free of every unknown. An empty result
    means the coefficient is constant in the unknowns.
    """
    unknowns = set(space.component_symbols)
    out = {}
    for sym in space.component_symbols:
        derivative = sympy.diff(expr, sym)
        if derivative != 0:
            out[sym.name] = "nonlinear" if derivative.free_symbols & unknowns else "linear"
    return out


def adr_form(equations: Sequence[sympy.Expr], space: Space,
             hamiltonian_gradients=()) -> Dict[str, Any]:
    """Classify scalar residuals and export them as plain data.

    The result is a dict that serializes to YAML. Per equation ``i`` (in
    component order) it holds, for each nonzero slot, the coefficient as a
    numpy code string over the names ``x, y, z, t``, the component names,
    and ``grad_<component>_<axis>`` (Hamiltonian only), plus its derivative
    with respect to each component it depends on, and ``depends_on``: how
    it depends on each (``linear`` or ``nonlinear``; an empty map is a
    constant). How a solver encodes that -- Proteus's flag dictionaries,
    say -- is the solver's business. Diffusion tensors are given in CSR
    form over their symbolic nonzeros.
    """
    if len(equations) != len(space.components):
        raise SymbolicError(
            "%d scalar equations for %d scalar unknowns (%s)"
            % (len(equations), len(space.components), ", ".join(space.components)))
    unknowns = set(space.component_symbols)
    gsubs = _gradient_symbols(space)
    comps = space.components
    out_equations = []
    for i, residual in enumerate(equations):
        eq = classify(residual, space, where="equation for %s" % comps[i],
                      hamiltonian_gradients=hamiltonian_gradients)
        entry: Dict[str, Any] = {"component": comps[i]}

        if eq.mass != 0:
            deps = _flag(eq.mass, space)
            if not deps:
                raise SymbolicError("equation for %s: the time derivative acts on no "
                                    "unknown" % comps[i])
            entry["mass"] = {
                "m": to_code(eq.mass),
                "dm": {c: to_code(sympy.diff(eq.mass, s)) for c, s in _syms(space, deps)},
                "depends_on": deps,
            }

        if any(f != 0 for f in eq.advection):
            deps = {}
            for f in eq.advection:
                for c, d in _flag(f, space).items():
                    deps[c] = "nonlinear" if "nonlinear" in (d, deps.get(c)) else d
            entry["advection"] = {
                "f": [to_code(f) for f in eq.advection],
                "df": {c: [to_code(sympy.diff(f, s)) for f in eq.advection]
                       for c, s in _syms(space, deps)},
                "depends_on": deps,
            }

        if eq.diffusion:
            diffusion = {}
            for phi, tensor in eq.diffusion.items():
                rows, cols, values = _csr(tensor)
                deps = {}
                for v in values:
                    for c, d in _flag(v, space).items():
                        deps[c] = "nonlinear" if "nonlinear" in (d, deps.get(c)) else d
                diffusion[phi] = {
                    "rowptr": rows, "colind": cols,
                    "a": [to_code(v) for v in values],
                    "da": {c: [to_code(sympy.diff(v, s)) for v in values]
                           for c, s in _syms(space, deps)},
                    "depends_on": deps,
                }
            entry["diffusion"] = diffusion

        if eq.reaction != 0:
            deps = _flag(eq.reaction, space)
            entry["reaction"] = {
                "r": to_code(eq.reaction),
                "dr": {c: to_code(sympy.diff(eq.reaction, s)) for c, s in _syms(space, deps)},
                "depends_on": deps,
            }

        if eq.hamiltonian != 0:
            h = eq.hamiltonian.xreplace(gsubs)
            dH = {}
            flags = {}
            for comp in comps:
                grads = [sympy.diff(h, sympy.Symbol(grad_name(comp, k), real=True))
                         for k in range(space.dim)]
                if any(g != 0 for g in grads):
                    dH[comp] = [to_code(g) for g in grads]
                    nonlinear = any(g.free_symbols & (set(gsubs.values()) | unknowns)
                                    for g in grads)
                    flags[comp] = "nonlinear" if nonlinear else "linear"
            if any(sympy.diff(h, s) != 0 for s in space.component_symbols):
                raise SymbolicError(
                    "equation for %s: the Hamiltonian %s depends on the unknowns "
                    "themselves, not only their gradients; Proteus's H(grad u) can't "
                    "hold that -- write the term in divergence form" % (comps[i], eq.hamiltonian))
            entry["hamiltonian"] = {"H": to_code(h), "dH": dH, "depends_on": flags}

        out_equations.append(entry)
    return {"dim": space.dim, "components": list(comps), "equations": out_equations}


def _syms(space: Space, deps: Dict[str, str]):
    return [(c, s) for c, s in zip(space.components, space.component_symbols) if c in deps]


def _csr(tensor: sympy.Matrix):
    rows, cols, values = [0], [], []
    for r in range(tensor.shape[0]):
        for c in range(tensor.shape[1]):
            v = sympy.simplify(tensor[r, c])
            if v != 0:
                cols.append(c)
                values.append(v)
        rows.append(len(cols))
    return rows, cols, values
