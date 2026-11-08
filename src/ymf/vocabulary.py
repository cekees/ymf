"""The words a spec states and a solver maps: ymf's vocabulary.

Each kind is a mapping ``{word: meaning}``. So far: :data:`LINEAR_SOLVERS`,
``direct`` and ``iterative`` and the packages a spec may still name, and
:func:`linear_solver_kind`. Stdlib only, like the rest of the spec core.
"""

from __future__ import annotations

from typing import Any, Dict, Optional

__all__ = ["LINEAR_SOLVERS", "linear_solver_kind", "VocabularyError"]


class VocabularyError(ValueError):
    """A word the vocabulary does not have, or one used where it cannot be."""


# ---------------------------------------------------------------------------
# Linear solvers
# ---------------------------------------------------------------------------

#: ``direct`` and ``iterative`` say what is asked of the linear solve; each
#: solver maps them to its own options and refuses what it cannot honour. A
#: package name says which implementation, and is a request of its kind.
LINEAR_SOLVERS: Dict[str, Dict[str, Any]] = {
    "direct": {"kind": "direct", "means": "a direct factorization; the solver picks the package "
                                          "(serial or distributed)"},
    "iterative": {"kind": "iterative", "means": "a Krylov method with a preconditioner; "
                                                "petsc_options refine it"},
    "superlu": {"kind": "direct", "means": "SuperLU (SuperLU_DIST in parallel)"},
    "mumps": {"kind": "direct", "means": "MUMPS"},
    "umfpack": {"kind": "direct", "means": "UMFPACK (serial)"},
    "petsc": {"kind": "iterative", "means": "a PETSc KSP; direct if petsc_options ask for "
                                            "pc_type lu or cholesky with ksp_type preonly"},
}

_DIRECT_PC = ("lu", "cholesky")


def linear_solver_kind(solver: Optional[Dict[str, Any]]) -> Optional[str]:
    """``direct``, ``iterative``, or None when the solver block names no
    linear solver."""
    if not solver or solver.get("linear_solver") is None:
        return None
    word = str(solver["linear_solver"])
    if word not in LINEAR_SOLVERS:
        raise VocabularyError("linear_solver %s: ymf knows %s" % (word, ", ".join(LINEAR_SOLVERS)))
    if word == "petsc":
        options = solver.get("petsc_options") or {}
        if (str(options.get("pc_type")) in _DIRECT_PC
                and str(options.get("ksp_type", "preonly")) == "preonly"):
            return "direct"
    return LINEAR_SOLVERS[word]["kind"]
