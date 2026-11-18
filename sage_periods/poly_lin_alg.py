r"""Linear algebra on polynomials and their encoded cofactors."""

from typing import NamedTuple

from sage.matrix.constructor import Matrix
from sage.matrix.special import identity_matrix

class ReducedElement(NamedTuple):
    r"""Pair a polynomial value with an optional encoded cofactor.

    A tracked cofactor is a polynomial or a ``CofactorCombination``; read it
    with ``as_polynomial``.

    INPUT:

    * ``value`` -- Encoded polynomial value.
    * ``cofactor`` -- (Optional) Encoded polynomial or ``CofactorCombination``;
      ``None`` means untracked.

    EXAMPLES:

        sage: from sage_periods.poly_lin_alg import ReducedElement
        sage: R = PolynomialRing(QQ, ['x', 'y']); x, y = R.gens()
        sage: ReducedElement(x).cofactor is None
        True
    """
    value: object
    cofactor: object | None = None


class CofactorCombination:
    r"""
    The cofactor $\sum_i c_i g_i$, stored as the coefficient vector $c$ over a
    list ``gens`` of input cofactors shared by one echelonization, and
    expanded only on demand.

    Expanding every echelon row would make each cofactor dense over all
    cofactor monomials of the slice; ``linear_normal_form_p`` instead combines
    the coefficient vectors and expands once per target.

    INPUT:

    * ``coeffs`` -- Coefficient vector.
    * ``gens`` -- List of encoded polynomial cofactors.

    EXAMPLES:

        sage: from sage_periods.poly_lin_alg import CofactorCombination
        sage: R = PolynomialRing(QQ, ['x', 'y']); x, y = R.gens()
        sage: c = CofactorCombination(vector(QQ, [2, 3]), [x, y])
        sage: c.polynomial()
        2*x + 3*y
    """
    __slots__ = ("coeffs", "gens", "_poly")

    def __init__(self, coeffs, gens):
        r"""Store ``coeffs`` over cofactor list ``gens``.

        INPUT:

        * ``coeffs`` -- Coefficient vector.
        * ``gens`` -- Encoded polynomial cofactors.

        TESTS:

            sage: from sage_periods.poly_lin_alg import CofactorCombination
            sage: c = CofactorCombination(vector(QQ, [1]), [PolynomialRing(QQ, ['x', 'y']).one()])
            sage: len(c.gens)
            1
        """
        self.coeffs, self.gens, self._poly = coeffs, gens, None

    def polynomial(self):
        r"""Expand combination to an encoded polynomial.

        OUTPUT:

        * The sum of ``coeffs[i] * gens[i]``.

        TESTS:

            sage: from sage_periods.poly_lin_alg import CofactorCombination
            sage: R = PolynomialRing(QQ, ['x', 'y']); x, y = R.gens()
            sage: CofactorCombination(vector(QQ, [1, 2]), [x, y]).polynomial()
            x + 2*y
        """
        if self._poly is None:
            self._poly = _combine(self.coeffs, self.gens)
        return self._poly


def _combine(coeffs, gens):
    r"""Return ``sum(coeffs[i] * gens[i])`` for a vector ``coeffs``."""
    out = gens[0].parent().zero()
    for i in coeffs.nonzero_positions():
        out += coeffs[i] * gens[i]
    return out


def as_polynomial(cofactor):
    r"""Return a tracked cofactor as a polynomial.

    INPUT:

    * ``cofactor`` -- An encoded polynomial or ``CofactorCombination``.

    OUTPUT:

    * The encoded polynomial.

    EXAMPLES:

        sage: from sage_periods.poly_lin_alg import CofactorCombination, as_polynomial
        sage: R = PolynomialRing(QQ, ['x', 'y']); x, y = R.gens()
        sage: as_polynomial(CofactorCombination(vector(QQ, [1, 2]), [x, y]))
        x + 2*y
    """
    if isinstance(cofactor, CofactorCombination):
        return cofactor.polynomial()
    return cofactor


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
    * ``return_pivot_rows`` -- (Optional) Boolean. If ``True``, also return
        the indices of input rows that pivot in the value block.

    ASSUMPTIONS: (not checked)

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

    EXAMPLES:

        sage: from sage_periods.poly_lin_alg import ReducedElement, echelonized_basis_poly
        sage: R = PolynomialRing(QQ, ["x", "y"])
        sage: x, y = R.gens()
        sage: L = [ReducedElement(x), ReducedElement(2*x)]
        sage: echelonized_basis_poly(L)
        [ReducedElement(value=x, cofactor=None)]
        sage: echelonized_basis_poly(L, return_matrix=True)
        [1]
        sage: echelonized_basis_poly(L, return_pivot_rows=True)[1]
        (0,)
    """
    # A ReducedElement with cofactor=None uses the same value-only matrix.
    if not L:
        return ([], []) if return_pivot_rows else []
    values = [item.value for item in L]
    # !! Assuming all cofactors are ReducedElements...
    cofactors = [item.cofactor for item in L]
    tracked = cofactors[0] is not None # Are we tracking cofactors or not?
    if tracked:
        # !! Assuming all cofactors are either all None or all Polys...
        cofactors = [as_polynomial(c) for c in cofactors]
    A = values[0].parent()

    # Enumerate all exponent vectors appearing in ``L`` and order them
    # decreasingly with respect to the monomial order on ``A``.  
    # It turns out working with exponent vecotrs directly is faster 
    # than with poly ring elts or a CombinatorialFreeModule.
    sortkey = A.term_order().sortkey
    dictionaries = [p.dict() for p in values]
    exps = set()
    for coefficients in dictionaries:
        exps.update(coefficients)
    exps = sorted(exps, key=sortkey, reverse=True)
    index = {e: j for j, e in enumerate(exps)}
    nvals = len(exps) # Number of value items

    if tracked and not return_matrix:
        # Echelonize [V | I]; the identity block records the row operations T.
        # If V has full row rank, T (with T*V = RREF(V)) is unique, so T*C is
        # the cofactor block of RREF([V | C]). Hence we keep T instead of T*C.
        F = A.base_ring()
        mat = Matrix(F, len(L), nvals, sparse=False)
        # Build V portion of matrix
        for i, coefficients in enumerate(dictionaries):
            for e, c in coefficients.items():
                mat[i, index[e]] = c
        piv = mat.pivot_rows() if return_pivot_rows else None
        # Augment identity matrix to transform cofactor columns.
        aug = mat.augment(identity_matrix(F, len(L)))
        aug.echelonize()
        
        # If we have only value pivots (i.e. V has full row rank), 
        # ... store [V | I] echelonization only; apply to CofactorCombination later, 
        # ... in linear_normal_form_p.
        if aug.pivots()[-1] < nvals: 
            reduced = aug.submatrix(0, 0, len(L), nvals)
            rows = [{} for _ in L]
            for i, j in reduced.nonzero_positions(copy=False):
                rows[i][exps[j]] = reduced[i, j]
            transform = aug.submatrix(0, nvals, len(L), len(L)).rows()
            out = [ReducedElement(A(row), CofactorCombination(coeffs, cofactors))
                   for row, coeffs in zip(rows, transform)]
            return (out, piv) if return_pivot_rows else out
        # Dependent values: use the [V | C] elimination below instead.

    # === Old version. Fall back to this when we don't have full row rank in V. ===

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

    # Rebuild matrix output
    # nonzero_positions() lists (row, column) pairs in row-major
    # order; group them by row instead of extracting each row vector.
    # Rows with a value pivot come first; stop at first cofactor-only row.
    out = []
    row_index, value, cofactor = None, None, None
    for i, j in mat.nonzero_positions(copy=False): # For each nonzero (row,col) pair (i,j)
        if i != row_index:  # Once we've gotten to a new row, we should append the ReducedElement from the previous row.
            if row_index is not None:
                out.append(ReducedElement(A(value), A(cofactor) if tracked else None))
            if j >= nvals: # "Recall" that nvals equals the number of value entries.
                row_index = None
                break  # A cofactor-only row supplies no reduction direction.
            row_index, value, cofactor = i, {}, {}
        if j < nvals:
            value[exps[j]] = mat[i, j]
        else:
            cofactor[cofactor_exps[j - nvals]] = mat[i, j]
    if row_index is not None: # Take care of the very last row after we exit the loop
        out.append(ReducedElement(A(value), A(cofactor) if tracked else None))
    return (out, piv) if return_pivot_rows else out


def linear_normal_form_p(L0, M, reduce):
    r"""
    Cancel an echelon basis's pivot monomials in each target.

    ``L0`` and ``M`` contain ``ReducedElement`` objects of the same cofactor
    mode. Value coefficients determine each subtraction and the same
    subtraction is applied to cofactors. If ``reduce`` is true,
    echelonize the residual family before returning it.

    When all of ``M`` carries ``CofactorCombination`` cofactors over one
    generator list, the subtractions are summed on the coefficient vectors
    and each target's cofactor is expanded once.

    INPUT:

    * ``L0`` -- A list of target ``ReducedElement`` objects, updated in place.
    * ``M`` -- An echelonized list of ``ReducedElement`` objects, representing an echelonized basis to reduce against.
    * ``reduce`` -- If ``True``, echelonize the residual list before returning it.

    OUTPUT:

    * The residual list of ``ReducedElement`` objects.

    EXAMPLES:

        sage: from sage_periods.poly_lin_alg import ReducedElement, linear_normal_form_p
        sage: R = PolynomialRing(QQ, ["x", "y"])
        sage: x, y = R.gens()
        sage: L = [ReducedElement(x + 1)]
        sage: linear_normal_form_p(L, [ReducedElement(x)], False)
        [ReducedElement(value=1, cofactor=None)]
        sage: L[0].value
        1
    """
    L = L0
    if not L:
        return []
    tracked = L[0].cofactor is not None
    # NOTE: Loop is formulated this way so we avoid building a new ReducedElement
    # after every single subtraction; only need one call per element of L0 at the end,
    # and this only if it is changed.
    gens = None
    if tracked and M and isinstance(M[0].cofactor, CofactorCombination):
        gens = M[0].cofactor.gens
        if any(not isinstance(m.cofactor, CofactorCombination)
               or m.cofactor.gens is not gens for m in M):
            gens = None # Sentinel allowing us to apply row op to entire CofactorCombination at once, or  proceed as before.
    pivots = [m.value.lm() for m in M]
    for i, p in enumerate(L):
        value, cofactor, changed, weights = p.value, p.cofactor, False, None
        for m, pivot in zip(M, pivots):
            c = value.monomial_coefficient(pivot)
            if c:
                value = value - c*m.value
                if gens is not None:
                    # Keep bulding batch transformation for cofactor row
                    w = c*m.cofactor.coeffs
                    weights = w if weights is None else weights + w
                elif tracked:
                    cofactor = as_polynomial(cofactor) - c*as_polynomial(m.cofactor)
                changed = True
        # Apply transformation to CofactorRow, if we have all cofactors stored in batch form
        if weights is not None:
            cofactor = as_polynomial(cofactor) - _combine(weights, gens)
        if changed:
            L[i] = ReducedElement(value, cofactor)
    
    # never entered
    if reduce:
        L = echelonized_basis_poly(L, False)

    return L
