r"""
Small linear algebra utilities for free modules over an algebra.

These routines are intended to work with Sage's
``CombinatorialFreeModule`` class and related structures.
"""
from sage.matrix.constructor import Matrix


# Future: added option to return pivot rows -- (possibly useful for profile)
def echelonized_basis_poly(L, return_matrix=False, return_pivot_rows=False):
    r"""
    Return a (reduced) echelonized basis for the span of a list of polynomials.
    That is, return a basis $B$ for span$_{\mathbb{K}}(L)$ such that for each $b \in B$ the leading monomial 
    of $b$ doesn't appear with nonzero coefficient in any other leading monomials of elements of $B$.

    This is an internal function for sage_periods, and is not 
    meant to be called by the user.

    INPUT:

    * ``L`` -- A list of elements of a multivariate polynomial ring ``A`` over a field ``K``.
    * ``return_matrix`` -- (Optional) A boolean with default value ``False``. 

    ASSUMPTIONS:

    * All elements in ``L`` have the same parent ``A``, a polynomial ring.

    OUTPUT:

    * A list of polynomials forming a basis for the $\mathbb{K}$-span of $L$ in
      echelon form. If ``return_matrix`` is ``True`` then a pair
      ``(basis, M)`` is returned, where ``M`` is the echelonized matrix
      whose rows correspond to the basis elements.

    NOTES:

    * ``L`` should consist of a list of vectors representing elements
        in a finite-dimensional ``A``-module.
    * The polynomial ring ``A`` should have grevlex order, which is
        compatible with the TOP order on our ``A``-module extending ``grevlex``.

    EXAMPLES:

    A nontrivial example.

            sage: var('t x0 x1 x2')
            sage: p = 7
            sage: A = PolynomialRing(GF(p),[x0,x1,x2])
            sage: x = A.gens()
            sage: f_arr = [ A.random_element() for i in range(5)]
            sage: f_arr
            [3*x0^2 + 2*x2^2 + 3*x0 - 2, -2*x0*x1 - 3*x1^2 - 3*x0*x2 - x2^2 - 3*x2, x0^2 - 3*x0*x1 - 3*x2^2 - x0, 3*x1^2 - 3*x0*x2 - 3*x2^2 + x0 - x1, -x0^2 - x0*x1 - x1*x2 + 3*x2^2 + x0]
            sage: echelonized_basis_poly(f_arr,False)
            [x0^2 + 3*x2^2 + x0 - 3, x0*x1 + 2*x2^2 + 3*x0 - 1, x1^2 - x2^2 - 2*x0 + x1 - 3*x2 - 2, x0*x2 - x1 - 3*x2 - 2, x1*x2 - x2^2 + 2*x0 - 3]
                                                                                                                                                                                                                                        0                                                                                                                                                                                                                                   0                                                                                                                                                                                                                                   0                                                                                                                                                                                                                                   0                                                                                                                                                                                                                                   1                                                                                                                                                                                            (-3/224*t^2 - 57/14)/(t^2 + 3/7*t + 3/7)                                                                                                                                                                                                           3/7*t/(t^2 + 3/7*t + 3/7)                                                                                                                                                                                                   (3/7*t + 6/7)/(t^2 + 3/7*t + 3/7)                                                                                                                                                                                                 (-9/14*t + 3/7)/(t^2 + 3/7*t + 3/7)                                                                                                                                                                                                                                   0]

    A trivial example.

            sage: echelonized_basis_poly([],False)
            []
            sage: echelonized_basis_poly([A(0)],False)
            []

    !!! note

            This method also works over $Q$, $Q(t)$, etc.
        
    """
    # Corner case / get A:
    if not L:
        return ([], []) if return_pivot_rows else []
    A = L[0].parent()

    # Enumerate all exponent vectors appearing in ``L`` and order them
    # decreasingly with respect to the monomial order on ``A``.  Working with
    # exponent vectors directly (rather than with ring elements or a
    # CombinatorialFreeModule) lets us build the coefficient matrix straight
    # from the polynomial dictionaries.
    sortkey = A.term_order().sortkey
    exps = set()
    for p in L:
        exps.update(p.dict().keys())
    exps = sorted(exps, key=sortkey, reverse=True)
    index = {e: j for j, e in enumerate(exps)}

    # Build the coefficient matrix and take RREF; this gives an echelonized
    # basis in coordinates.  Over GF(p) this uses fast native linear algebra.
    entries = {}
    for i, p in enumerate(L):
        for e, c in p.dict().items():
            entries[(i, index[e])] = c
    mat = Matrix(A.base_ring(), len(L), len(exps), entries, sparse=False)

    # If requested, record which input rows form a basis of the row space
    # (this must be read off before echelonizing, which mutates the matrix).
    piv = mat.pivot_rows() if return_pivot_rows else None

    mat.echelonize()

    if return_matrix:
        rows = [row for row in mat.rows() if not row.is_zero()]
        return (Matrix(rows), piv) if return_pivot_rows else Matrix(rows)

    # Convert the nonzero rows back to polynomials.  RREF places zero rows
    # last, so we may stop at the first one.
    out = []
    for row in mat.rows():
        nz = row.nonzero_positions()
        if not nz:
            break
        out.append(A({exps[j]: row[j] for j in nz}))
    return (out, piv) if return_pivot_rows else out
    

def linear_normal_form_p(L0 , M, reduce):
    r"""
    Reduce multivariate polynomials modulo an already-echelonized list, directly in polynomial form.
    
    We assume $M$ is already echelonized, so reduction is performed simply by cancelling each leading 
    monomial of $M$ in all entries of $L_0$.

    This is an internal function for sage_periods, and is not 
    meant to be called by the user.
    
    INPUT:
    
    * ``L0`` -- A sequence of multivariate polynomials to reduce.
    * ``M`` -- An echelonized sequence of multivariate polynomials.
    * ``reduce`` -- A boolean flag signalling whether to echelonize the output again.
    
    OUTPUT:
    
    * A sequence of polynomial normal forms modulo the span of M.

    EXAMPLES:

    A nontrivial example.

            sage: var('t x0 x1 x2')
            sage: p = 7
            sage: A = PolynomialRing(GF(p),[x0,x1,x2])
            sage: x = A.gens()
            sage: B = echelonized_basis_poly([ A.random_element() for i in range(3)],A,False)
            sage: linear_normal_form_p([A.random_element() for i in range(2)],B,reduce=False)
            [-x0^2 - x1 + x2 - 3, -2*x0^2 + 2*x2^2 - 2*x0 + 2*x2 + 3]

    You can reduce the result with respect to ``M0`` if you want to.

            sage: linear_normal_form_p([A.random_element() for i in range(2)],B,reduce=True)
            [x0*x1 - x1 + 3*x2 - 2, x2^2 - 3*x2 - 3]

    Trivial examples.

            sage: linear_normal_form_p([],B,False) # Trivial input list
            []
            sage: linear_normal_form_p([A(3*x[2]+1)],[],False) # Trivial basis reduction gives itself.
            [3*x2 + 1]

    """
    L = L0
    # Successively cancel the pivont monomial of each row of M in all targets
    for j in range(len(M)):
        m = M[j].lm()
        for i in range(len(L)):
            c = L[i].monomial_coefficient(m)
            L[i] -= c*M[j]

    if reduce:
        # Optinally normalize the residual family itself
        L = echelonized_basis_poly(L,False)
    
    return L