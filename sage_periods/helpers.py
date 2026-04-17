r"""
Holds basic helper functions for the package.
"""

from sage.rings.number_field.number_field import CyclotomicField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.symbolic.ring import SR
from sage.misc.misc_c import prod


def _compress_invariant_poly(P, idx, m, Anew):
    """
    P is a polynomial invariant under x_idx -> zeta * x_idx.
    Therefore every exponent of variable idx is divisible by $m$.
    Replace x_idx^(m*k) by newvar^k.
    """
    new_gens = Anew.gens()
    out = Anew.zero()

    # Can't just use subs; best to deal with exponent vectors.
    for exp, c in P.dict().items():
        e = exp[idx]
        if e % m != 0:
            raise ValueError("Polynomial is not m-invariant in the chosen variable.")
        new_exp = list(exp)
        new_exp[idx] = e // m
        mon = prod(new_gens[i]**new_exp[i] for i in range(len(new_exp)))
        out += Anew(c) * mon

    return out


def _coerce_poly_to_QQ(P, AQQ):
    """
    Try to coerce all coefficients of P exactly into QQ.
    """
    qq_gens = AQQ.gens()
    out = AQQ.zero()

    for exp, c in P.dict().items():
        cq = QQ(c)   # raises if coefficient is not rational
        mon = prod(qq_gens[i]**exp[i] for i in range(len(exp)))
        out += cq * mon

    return out


def first_coordinate_section(R, z1, m, u=None, vari=None):
    r"""
    Takes a symbolic rational function and returns its first-coordinate
    section, which is the subseries including only first-index entries
    which are divisible by $m$. This is obtained by applying the $m$th root
    of unity filter in the first coordinate to $R$.

    INPUT:

    * ``R``   -- A symbolic rational function.
    * ``z1``  -- The variable to section.
    * ``m``   -- A positive integer.
    * ``u``    -- (Optional) A replacement variable after sending $z_1^{mk} \mapsto u^k$. Defaults to ``z1``.
    * ``vari`` -- (Optional) An ordered variable list. Defaults to ``sorted(R.variables())``.

    ASSUMPTIONS:

    * ``R`` has coefficients in ``QQ`` (or at least in a field that survives coercion through the chosen cyclotomic field).
    * ``m`` is greater than 1. If ``m`` is 0 then our code in ``picard_fuchs`` uses a different formula for computing its diagonal.

    OUTPUT:

    * ``out`` -- A rational function equal to the first-coordinate ``m``-section of ``R``.

    EXAMPLES:

    A nontrivial example.

            sage: var('t x y')
            sage: F = 1/(1-x-y-x*y**3)
            sage: r=[2,3]
            sage: first_coordinate_section(F,x, 2)
            (y - 1)/(x*y^6 + 2*x*y^3 - y^2 + x + 2*y - 1)

    A trivial example.

            sage: # If a constant is passed in, an error will be thrown.
            sage: first_coordinate_section(SR(1),x, 7)
            Traceback (most recent call last):
            ...
            ValueError: z1 must be one of the variables of R.

    """

    if u is None:
        u = z1
    if vari is None:
        vari = sorted(R.variables(), key=str)
    if z1 not in vari:
        raise ValueError("z1 must be one of the variables of R.")
    
    # Assuming m > 1
    names_old = [str(v) for v in vari]
    idx = names_old.index(str(z1))

    # Work over a cyclotomic field containing a primitive m-th root of unity
    Kz = CyclotomicField(m)
    zeta = Kz.gen()

    A = PolynomialRing(Kz, names_old, order="degrevlex")
    B = A.fraction_field()
    gens = A.gens()
    xx = gens[idx]

    # Convert R into B
    RB = B(R)
    N = RB.numerator()
    D = RB.denominator()

    # Build the roots-of-unity average with a common denominator:
    #   (1/m) sum_j N(zeta^j x)/D(zeta^j x)
    # = (1/m) * Num / Den

    # This common denominator is necessary
    Njs = []
    Djs = []
    for j in range(m):
        subs_dict = {gens[k]: gens[k] for k in range(len(gens))}
        subs_dict[xx] = (zeta**j) * xx
        Njs.append(A(N.subs(subs_dict)))
        Djs.append(A(D.subs(subs_dict)))

    Den = prod(Djs)
    Num = A.zero()
    for j in range(m):
        term = Njs[j]
        for i in range(m):
            if i != j:
                term *= Djs[i]
        Num += term
    Num *= Kz(1) / Kz(m)

    # Rewrite z1^(m*k) -> u^k
    names_new = [str(u) if k == idx else names_old[k] for k in range(len(names_old))]
    Anew = PolynomialRing(Kz, names_new, order="degrevlex")

    Num_c = _compress_invariant_poly(Num, idx, m, Anew)
    Den_c = _compress_invariant_poly(Den, idx, m, Anew)

    # Put coefficients back into QQ (We should always be able to do this,
    # because first-coordinate section always has rational coefficients.)
    AQQ = PolynomialRing(QQ, names_new, order="degrevlex")
    Num_q = _coerce_poly_to_QQ(Num_c, AQQ)
    Den_q = _coerce_poly_to_QQ(Den_c, AQQ)
    out = AQQ.fraction_field()(Num_q) / AQQ.fraction_field()(Den_q)

    return SR(out)