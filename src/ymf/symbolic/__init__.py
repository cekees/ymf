"""The symbolic layer: machine-readable equations -> solver-ready forms.

Needs sympy (``pip install 'ymf[symbolic]'``). Its output is plain data, so
solvers that consume it do not need sympy; see :mod:`ymf.symbolic.adr`.
"""

from ymf.symbolic.adr import adr_form, classify, grad_name, to_code
from ymf.symbolic.language import (D, Dt, Space, SymbolicError, parse_equation,
                                   parse_expression)
from ymf.symbolic.problem import adr_problem

__all__ = ["D", "Dt", "Space", "SymbolicError", "parse_equation", "parse_expression",
           "classify", "adr_form", "grad_name", "to_code", "adr_problem"]
