"""linear_solver: direct / iterative (C42 round 3, decision 19): what is
asked of the linear solve, shared by every solver path; a package name is a
request of its kind.
"""

import pytest

from ymf.vocabulary import VocabularyError, linear_solver_kind


@pytest.mark.parametrize("solver, kind", [
    ({"linear_solver": "direct"}, "direct"),
    ({"linear_solver": "iterative"}, "iterative"),
    ({"linear_solver": "superlu"}, "direct"),
    ({"linear_solver": "mumps"}, "direct"),
    ({"linear_solver": "umfpack"}, "direct"),
    ({"linear_solver": "petsc"}, "iterative"),
    ({"linear_solver": "petsc", "petsc_options": {"ksp_type": "gmres", "pc_type": "asm"}}, "iterative"),
    ({"linear_solver": "petsc", "petsc_options": {"ksp_type": "preonly", "pc_type": "lu"}}, "direct"),
    ({"linear_solver": "petsc", "petsc_options": {"pc_type": "cholesky"}}, "direct"),
    ({"linear_solver": "petsc", "petsc_options": {"ksp_type": "cg", "pc_type": "lu"}}, "iterative"),
    ({"nonlinear_solver": "newton"}, None),
    (None, None),
])
def test_the_kind_of_a_linear_solve(solver, kind):
    assert linear_solver_kind(solver) == kind


def test_an_unknown_linear_solver_is_refused():
    with pytest.raises(VocabularyError, match="linear_solver hypre: ymf knows direct, iterative"):
        linear_solver_kind({"linear_solver": "hypre"})

