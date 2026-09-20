# cython: language_level=3, overflowcheck=True

from sage.rings.rational_field import QQ
from sage.rings.rational cimport Rational


def diagonal_series_terms_cython(P, Q, r, Py_ssize_t N):
    r"""
    Return the first ``N`` terms of the ``r``-diagonal of the power series
    expansion of $P/Q$ at the origin, where $P$ and $Q$ are polynomials over
    $\mathbb{Q}$ in the same variables and $Q(0) \neq 0$.

    The expansion of $1/Q$ is computed layer by layer in the total degree,
    pruning every exponent outside the box $\prod_i [0, (N-1) r_i]$: since
    polynomial multiplication only increases exponents, no pruned monomial can
    contribute to a diagonal coefficient inside the box, so the computed terms
    are exact.

    Loop counters and coordinate arithmetic use checked machine-sized integers;
    coefficients remain arbitrary-precision Sage rationals.

    This is an internal function for sage_periods, and is not meant to be
    called by the user.
    """
    cdef Py_ssize_t d = P.parent().ngens()
    cdef Py_ssize_t m, j, i, de, maxdeg, a, b, coordinate
    cdef tuple zero, box, e, ey, enew, target, rem
    cdef dict Qd, Pd, layers, layer, src, ycoeff, lay
    cdef list Qterms, terms, coords
    cdef Rational c0, c, cy, s, v
    cdef Rational qzero = QQ(0)

    assert len(r) == d, "r must have same length as number of generators of P and Q."
    zero = (0,) * d
    Qd = {tuple(raw_e): QQ(raw_c) for raw_e, raw_c in Q.dict().items()}
    c0 = Qd.pop(zero, None)
    if c0 is None or c0 == 0:
        raise ValueError("The denominator vanishes at the origin, so the power series diagonal is not defined.")
    box = tuple((N - 1) * ri for ri in r)
    maxdeg = sum(box)
    coords = [0] * d

    # Homogeneous-layer recurrence for y = 1/Q: writing Q = c0 + (higher order),
    # the degree-m part of y is y_m = -(1/c0) * sum_{|e|>=1} Q_e * y_{m-|e|}.
    layers = {0: {zero: 1 / c0}}
    Qterms = []
    for e, c in Qd.items():
        for i in range(d):
            if e[i] > box[i]:
                break
        else:
            Qterms.append((e, -c / c0, sum(e)))

    for m in range(1, maxdeg + 1):
        layer = {}
        for e, c, de in Qterms:
            src = layers.get(m - de)
            if not src:
                continue
            for ey, cy in src.items():
                for i in range(d):
                    a = ey[i]
                    b = e[i]
                    coordinate = a + b
                    if coordinate > box[i]:
                        break
                    coords[i] = coordinate
                else:
                    enew = tuple(coords)
                    layer[enew] = layer.get(enew, qzero) + c * cy
        layer = {e: c for e, c in layer.items() if c}
        if layer:
            layers[m] = layer
    ycoeff = {}
    for lay in layers.values():
        ycoeff.update(lay)

    # Multiply by P and extract the diagonal coefficients;
    # standard Cauchy product.
    Pd = {tuple(raw_e): QQ(raw_c) for raw_e, raw_c in P.dict().items()}
    terms = []
    for j in range(N):
        target = tuple(j * ri for ri in r)
        s = qzero
        for e, c in Pd.items():
            for i in range(d):
                a = target[i]
                b = e[i]
                coordinate = a - b
                if coordinate < 0:
                    break
                coords[i] = coordinate
            else:
                rem = tuple(coords)
                v = ycoeff.get(rem)
                if v is not None:
                    s += c * v
        terms.append(s)
    return terms
