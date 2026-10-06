"""The equation language: unevaluated operators, and parsing into sympy.

Equations are written in a small operator syntax that is ordinary Python
expression syntax, evaluated by sympy's parser in a namespace that holds
only the operators, the elementary functions, the coordinates, the unknowns
and the coefficients::

    rho*dt(v) + div(rho*outer(v, v)) - div(mu*grad(v)) + grad(p) = f
    div(v) = 0

The design follows ibvp (https://github.com/ibvp/ibvp, MIT, by Andreas
Kloeckner and Robert C. Kirby): fields, the operators dt, div, grad, and a
*scalarization* step that writes every vector equation as one scalar
equation per component, in terms of a single first-derivative operator.
Here that operator is :class:`D`, and the time derivative is :class:`Dt`.
Both stay unevaluated when applied to an unknown -- that structure is what
:mod:`ymf.symbolic.adr` reads to recover the flux form of each term -- and
evaluate to ordinary derivatives when applied to a known function of x and
t, so source terms and coefficients differentiate as expected.

A vector is a sympy column ``Matrix`` with one row per space dimension,
written ``(a, b)`` anywhere in an expression; a rank-2 tensor (``grad`` of
a vector, ``outer``) is a square ``Matrix``. Tensors follow the continuum
mechanics convention ``grad(v)[i, j] = d v_j / d x_i``, and ``div``
contracts the first index, so ``dot(v, grad(v))`` is (v . nabla) v and
``div(mu*grad(v))`` is the vector Laplacian.
"""

from __future__ import annotations

from typing import Dict, Iterable, List, Sequence

import sympy

__all__ = [
    "D", "Dt", "Space", "SymbolicError", "parse_expression", "parse_equation",
]


class SymbolicError(ValueError):
    """A machine-readable expression in a spec can't be understood."""


class D(sympy.Function):
    """The derivative of an unknown along one axis: ``D(u, 1)`` is du/dy.

    Constructed through :meth:`Space.d`, which keeps D linear and applies it
    only to the unknowns; known expressions are differentiated outright.
    """
    nargs = 2


class Dt(sympy.Function):
    """The time derivative of an expression in the unknowns."""
    nargs = 1


class Space:
    """The coordinates, unknowns and coefficients one problem is written in.

    ``unknowns`` maps a name to its rank: 0 for a scalar, 1 for a vector.
    ``coefficients`` maps a name to a sympy value (a number, or an
    expression in the coordinates and time); they are substituted as the
    equations are parsed, so a derivative of a variable coefficient comes
    out right.
    """

    def __init__(self, dim: int, unknowns: Dict[str, int],
                 coefficients: Dict[str, sympy.Expr] | None = None):
        if dim not in (1, 2, 3):
            raise SymbolicError("space dimension must be 1, 2 or 3, not %r" % (dim,))
        self.dim = dim
        self.x = sympy.symbols("x y z", real=True)[:dim]
        self.t = sympy.Symbol("t", real=True)
        self.unknowns = dict(unknowns)
        #: name of each scalar component, in order: "u", or "v_0", "v_1", ...
        self.components: List[str] = []
        self._fields: Dict[str, object] = {}
        for name, rank in self.unknowns.items():
            if rank == 0:
                self.components.append(name)
                self._fields[name] = sympy.Symbol(name, real=True)
            elif rank == 1:
                names = ["%s_%d" % (name, i) for i in range(dim)]
                self.components.extend(names)
                self._fields[name] = sympy.Matrix([sympy.Symbol(n, real=True) for n in names])
            else:
                raise SymbolicError("unknown %r has rank %r; only 0 and 1 are supported"
                                    % (name, rank))
        self.component_symbols = [sympy.Symbol(n, real=True) for n in self.components]
        self.coefficients = dict(coefficients or {})

    # -- operators ---------------------------------------------------------

    def is_known(self, expr) -> bool:
        return not (sympy.sympify(expr).free_symbols & set(self.component_symbols))

    def d(self, expr, axis: int):
        """d(expr)/dx_axis, linear, unevaluated only where unknowns appear."""
        expr = sympy.sympify(expr)
        if self.is_known(expr) and not expr.has(D, Dt):
            return sympy.diff(expr, self.x[axis])
        if isinstance(expr, sympy.Add):
            return sympy.Add(*[self.d(a, axis) for a in expr.args])
        if isinstance(expr, sympy.Mul):
            numbers = [f for f in expr.args if f.is_number]
            if numbers:
                rest = sympy.Mul(*[f for f in expr.args if not f.is_number])
                return sympy.Mul(*numbers) * self.d(rest, axis)
        return D(expr, axis)

    def dt(self, expr):
        if isinstance(expr, sympy.MatrixBase):
            return expr.applyfunc(self.dt)
        expr = sympy.sympify(expr)
        if self.is_known(expr):
            return sympy.diff(expr, self.t)
        if isinstance(expr, sympy.Add):
            return sympy.Add(*[self.dt(a) for a in expr.args])
        if isinstance(expr, sympy.Mul):
            numbers = [f for f in expr.args if f.is_number]
            if numbers:
                rest = sympy.Mul(*[f for f in expr.args if not f.is_number])
                return sympy.Mul(*numbers) * self.dt(rest)
        return Dt(expr)

    def grad(self, expr):
        if isinstance(expr, sympy.MatrixBase):
            if expr.shape[1] != 1:
                raise SymbolicError("grad of a %dx%d tensor is not supported" % expr.shape)
            return sympy.Matrix(self.dim, self.dim,
                                lambda i, j: self.d(expr[j], i))
        return sympy.Matrix([self.d(expr, j) for j in range(self.dim)])

    def div(self, expr):
        if not isinstance(expr, sympy.MatrixBase):
            raise SymbolicError("div needs a vector or a tensor, got the scalar %s" % (expr,))
        if expr.shape == (self.dim, 1):
            return sympy.Add(*[self.d(expr[j], j) for j in range(self.dim)])
        if expr.shape == (self.dim, self.dim):
            return sympy.Matrix([sympy.Add(*[self.d(expr[i, j], i) for i in range(self.dim)])
                                 for j in range(self.dim)])
        raise SymbolicError("div of a %dx%d object in %d dimensions" % (expr.shape + (self.dim,)))

    def dot(self, a, b):
        """Contract the last index of a with the first of b."""
        if isinstance(a, sympy.MatrixBase) and isinstance(b, sympy.MatrixBase):
            a_vec, b_vec = a.shape[1] == 1, b.shape[1] == 1
            if a_vec and b_vec:
                return (a.T * b)[0, 0]
            if a_vec:
                return (a.T * b).T
            return a * b
        raise SymbolicError("dot needs two vectors, or a vector and a tensor")

    def outer(self, a, b):
        if not (isinstance(a, sympy.MatrixBase) and isinstance(b, sympy.MatrixBase)):
            raise SymbolicError("outer needs two vectors")
        return a * b.T

    def namespace(self) -> Dict[str, object]:
        """Everything an equation string may name."""
        ns: Dict[str, object] = {
            "dt": self.dt, "div": self.div, "grad": self.grad,
            "dot": self.dot, "outer": self.outer,
            "lap": lambda e: self.div(self.grad(e)),
            "identity": sympy.eye(self.dim),
            "pi": sympy.pi, "E": sympy.E,
            "sin": sympy.sin, "cos": sympy.cos, "tan": sympy.tan,
            "exp": sympy.exp, "log": sympy.log, "sqrt": sympy.sqrt,
            "sinh": sympy.sinh, "cosh": sympy.cosh, "tanh": sympy.tanh,
            "abs": sympy.Abs, "Abs": sympy.Abs,
            "vec": lambda *xs: sympy.Matrix([sympy.sympify(v) for v in xs]),
            "t": self.t,
        }
        for axis, name in enumerate("xyz"[:self.dim]):
            ns[name] = self.x[axis]
            ns["d" + name] = (lambda e, a=axis: self.d(e, a))
        ns.update(self.coefficients)
        ns.update(self._fields)
        return ns


def parse_expression(text: str, space: Space):
    """Parse one expression in the spec notation (see :mod:`ymf.symbolic.notation`)."""
    from ymf.symbolic.notation import NotationError, Parser
    try:
        return Parser(space).expression(str(text))
    except NotationError as exc:
        raise SymbolicError(str(exc)) from exc


def parse_equation(text: str, space: Space) -> List[sympy.Expr]:
    """Parse ``lhs = rhs [in region]`` into scalar residuals ``lhs - rhs``."""
    from ymf.symbolic.notation import NotationError, Parser
    try:
        return Parser(space).equation(str(text))
    except NotationError as exc:
        raise SymbolicError(str(exc)) from exc
