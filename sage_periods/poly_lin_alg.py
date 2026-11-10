r"""Linear algebra on polynomials and their encoded cofactors."""

from typing import NamedTuple

from sage.matrix.constructor import Matrix

class ReducedElement(NamedTuple):
    """Value and cofactor in ``U.xring``; ``None`` means untracked."""
    value: object
    cofactor: object | None = None


def echelonized_basis_poly(L, return_matrix=False, return_pivot_rows=False):
    r"""
    Return a (reduced) echelonized basis for a list of ``ReducedElement`` objects.
    That is:
    * For ``ReducedElement`` objects, use only the ``value`` columns as pivots,
      while applying row operations to both components.

    INPUT:

    * ``L`` -- A list of ``ReducedElement`` objects.
    * ``return_matrix`` -- (Optional) Boolean. If ``True``, return their 
        coefficient matrix (value columns first) instead of the element list.
    * ``return_pivot_rows`` -- (Optional) Boolean If ``True``, Return 
        the input rows that pivot in the value block.

    ASSUMPTIONS:

    * All values are in one polynomial ring over a field. In the modular
      reduction pipeline this is ``U.xring = GF(p)[u,v,x_0,...]``.
    * A tracked cofactor is an encoded n-form in that same ring, not a list
      of coefficients in the original polynomial ring.
    * Inputs are all ``ReducedElement`` objects.

    OUTPUT:

    * A list in the same format as the input, with values in echelon form.
      With ``return_matrix=True``, return the coefficient matrix instead.
      With ``return_pivot_rows=True``, also return the input pivot rows.

    This is an internal function for sage_periods, and is not 
    meant to be called by the user.
    """
    # A ReducedElement with cofactor=None uses the same value-only matrix.
    if not L:
        return ([], []) if return_pivot_rows else []
    if any(not isinstance(item, ReducedElement) for item in L):
        raise TypeError("All inputs must be ReducedElement objects.")
    values = [item.value for item in L]
    cofactors = [item.cofactor for item in L]
    tracked = cofactors[0] is not None
    if any((c is not None) != tracked for c in cofactors):
        raise ValueError("Cofactors must either all be present or all be None.")
    A = values[0].parent()
    if any(p.parent() != A for p in values) or \
            (tracked and any(c.parent() != A for c in cofactors)):
        raise TypeError("Values and cofactors must belong to the same polynomial ring.")

    # Enumerate all exponent vectors appearing in ``L`` and order them
    # decreasingly with respect to the monomial order on ``A``.  Working with
    # exponent vectors directly (rather than with ring elements or a
    # CombinatorialFreeModule) lets us build the coefficient matrix straight
    # from the polynomial dictionaries.
    sortkey = A.term_order().sortkey
    dictionaries = [p.dict() for p in values]
    exps = set()
    for coefficients in dictionaries:
        exps.update(coefficients)
    exps = sorted(exps, key=sortkey, reverse=True)
    index = {e: j for j, e in enumerate(exps)}

    cofactor_dicts = [c.dict() for c in cofactors] if tracked else []
    cofactor_exps = sorted({e for d in cofactor_dicts for e in d},
                           key=sortkey, reverse=True)
    cofactor_index = {e: len(exps) + j for j, e in enumerate(cofactor_exps)}

    # Build the coefficient matrix and take RREF; this gives an echelonized
    # basis in coordinates.  Over GF(p) this uses fast native linear algebra.
    # Fill a zero dense matrix entry by entry: much faster than passing a dict
    # to the Matrix constructor.
    mat = Matrix(A.base_ring(), len(L), len(exps) + len(cofactor_exps), sparse=False)
    for i, coefficients in enumerate(dictionaries):
        for e, c in coefficients.items():
            mat[i, index[e]] = c
    for i, coefficients in enumerate(cofactor_dicts):
        for e, c in coefficients.items():
            mat[i, cofactor_index[e]] = c

    # Profiling only cares about the value rows, not cofactor-only rows.
    piv = None
    if return_pivot_rows:
        piv = mat[:, :len(exps)].pivot_rows() if tracked else mat.pivot_rows()

    mat.echelonize()

    # return_matrix==True case:
    if return_matrix:
        rows = []
        for row in mat.rows():
            nz = row.nonzero_positions()
            if not nz or nz[0] >= len(exps):
                break
            rows.append(row)
        return (Matrix(rows), piv) if return_pivot_rows else Matrix(rows)

    out = []
    for row in mat.rows():
        nz = row.nonzero_positions()
        if not nz or nz[0] >= len(exps):
            break  # A cofactor-only row supplies no reduction direction.
        value = A({exps[j]: row[j] for j in nz if j < len(exps)})
        cofactor = (A({cofactor_exps[j - len(exps)]: row[j]
                       for j in nz if j >= len(exps)}) if tracked else None)
        out.append(ReducedElement(value, cofactor))
    return (out, piv) if return_pivot_rows else out


def linear_normal_form_p(L0, M, reduce):
    r"""
    Cancel an echelon basis's pivot monomials in each target.

    ``L0`` and ``M`` contain ``ReducedElement`` objects of the same cofactor
    mode. Value coefficients determine each subtraction and the same
    subtraction is applied to cofactors. If ``reduce`` is true,
    echelonize the residual family before returning it.
    """
    L = L0
    if not L:
        return []
    tracked = L[0].cofactor is not None
    for m in M:
        pivot = m.value.lm()
        for i, p in enumerate(L):
            c = p.value.monomial_coefficient(pivot)
            if c:
                cofactor = p.cofactor - c*m.cofactor if tracked else None
                L[i] = ReducedElement(p.value - c*m.value, cofactor)

    if reduce:
        L = echelonized_basis_poly(L, False)

    return L
