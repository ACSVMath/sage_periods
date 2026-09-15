r"""Utilities for working with diagonals."""

from sage.all import *

def normalize_diagonal_arguments(R, r = None, vari = None, Dt = None, t = None):
    r"""
        Runs input validation of arguments passed into compute_diagonal_annihilator.
        Then it makes sure that all of R,r,vari,Dt,t are defined and of the correct format to be processed by diagonal_to_period.

        TODO
    """

    # Basics checks and argument processing
    # assert Dt != 'D', "Please pick a different operator symbol; D is reserved."
    assert Dt == None or isinstance(Dt,str), "Derivative symbol Dt must be None or a string."
    assert not(Dt != None and t == None), "If you pass Dt, then you must pass t as well."
    if t != None and Dt == None:
        Dt = f"D{t}"

    # Build "R.variables()" surrogate if we're not in SR
    if R.parent() != SR:
        Rvariables_poly = sum(set(R.numerator().variables()).union(R.denominator().variables())).monomials()

    # Default behavior if the user does not pass t nor Dt
    if t == None and Dt == None:
        t = SR.var('t')
        Dt = 'Dt'

    # Default behavior if the user does not pass r
    if r != None:
        d = len(r)
        if vari == None:
            if R.parent() == SR:
                print(f"WARNING: You specified a direction vector but not a list of variables. The ordering {R.variables()} will be used")
            else:
                print(f"WARNING: You specified a direction vector but not a list of variables. The ordering {Rvariables_poly} will be used")
        else:
            assert len(r) == len(vari), "Direction vector r must have same length as the number of variables."
        assert all(isinstance(x,int) or isinstance(x,Integer) for x in r), "Direction vector r must be a list of integers."
        assert (0 not in r) and (Integer(0) not in r), "Cannot have zero entry in r; in any case, this is isomorphic to the case in d-1 variables."
        assert all( ri > 0 or ri > Integer(0) for ri in r), "Cannot have a negative integer coordinate in r."
    else:
        if R.parent() == SR: 
            d = len(R.variables())
        else:
            d = len(Rvariables_poly)
        r = [1]*d

    if vari != None:
        assert all(x.parent() == SR for x in vari) or all(x in R.parent().gens() for x in vari), "Variable list must contain all symbolic variables or generators of R's parent."
        assert (set(vari) == set(R.variables())  or set(vari) == set(Rvariables_poly)) and len(vari) == d, "vari must contain exactly the variables appearing in R (except t)."
    else:
        if R.parent() == SR:
            vari = R.variables()
        else:
            vari = [ SR(v) for v in Rvariables_poly]
    
    # !! UNNECESSARY
    # # So as to not have to case on whether our input is symbolic or not when creating / destroying symbolic variables, just normalize everything to new symbolic variable names.
    # vari_normalized = list(SR.var("sage_periods_x", len(vari)))
    # R_normalized = SR(R)
    # R_normalized = R_normalized.subs({old:new for (old,new) in zip([ SR(v) for v in vari],vari_normalized)})

    R = SR(R)
    assert t not in R.variables()

    # Permute r, vari, vari_normalized so that smallest r-entry lies in front, saving computation in first_coordinate_section.
    _min_idx = r.index(min(r))
    _indices_permuted = [_min_idx] + [i for i in range(d) if i != _min_idx]
    r = [r[i] for i in _indices_permuted]
    vari = [vari[i] for i in _indices_permuted]

    return R, r, vari, Dt, t

def diagonal_to_period(R, r, vari , t):
    # Takes symbolic R and direction-variable order pair (r,vari) and converts this to the period integrand of the corresponding diagonal in t,
    # ...applying a RoU filter, if necessary.

    # !! Assumes R, r, vari, t are normalized and that the smallest entry of r lies in front.
    d = len(r)
    if r[0] > 1:
        # Build root of unity filter / first-coordinate section of F over cyclotomic field over $\zeta$,
        # where $\zeta$ is a primitive r1-th root of unity.
        m = r[0]
        Rsec = first_coordinate_section(R, vari[0], m, u=vari[0], vari=vari)
        chvar = [vari[0] / prod(vari[i]**(m * r[i]) for i in range(1, d))] \
                + [vari[i]**m for i in range(1, d)]
        G = Rsec.subs({vari[i]: chvar[i] for i in range(d)}) / prod(vari[1:])
    else:
        # r[0] = 1 means we can use the usual change of variables formula for G, and sub into F directly.
        chvar = [vari[0]/prod(vari[i]**(r[i]) for i in range(1,d))] + [vari[i] for i in range(1,d)]
        G = R.subs({vari[i]: chvar[i] for i in range(d)})/prod(vari[1:]) 
            
    # Change the first variable of G to t. If t is provided and is in F, but is not the first variable, then we swap t and vari[0].
    if t in vari and t != vari[0]:
        G = G.subs({ vari[0] : t, t:vari[0]})
    else:
        G = G.subs({ vari[0] : t})

    return G

def minimize_diagonal_annihilator(L, R, r, t):
    r"""
        TODO: write this function which will:
        - Compute initial terms in the diagonal of R = P/Q by series expansion.
        - "Guess" recurrence operator based on these initial terms.
        - Verify that this operator is correct, and is the minimal annihilator of the period.
    """
    pass

# Can get from AI version picard_fuchs.py, iyw.
def diagonal_series_terms(P, Q, rdir, N):
    r"""
    Return the first ``N`` terms of the ``rdir``-diagonal of the power series
    expansion of $P/Q$ at the origin, where $P$ and $Q$ are polynomials over
    $\mathbb{Q}$ in the same variables and $Q(0) \neq 0$.

    The expansion of $1/Q$ is computed layer by layer in the total degree,
    pruning every exponent outside the box $\prod_i [0, (N-1) r_i]$: since
    polynomial multiplication only increases exponents, no pruned monomial can
    contribute to a diagonal coefficient inside the box, so the computed terms
    are exact.

    This is an internal function for sage_periods, and is not meant to be
    called by the user.
    """
    d = P.parent().ngens()
    zero = tuple([0] * d)
    Qd = {tuple(e): QQ(c) for e, c in Q.dict().items()}
    c0 = Qd.pop(zero, None)
    if c0 is None or c0 == 0:
        raise ValueError("The denominator vanishes at the origin, so the power series diagonal is not defined.")
    box = tuple((N - 1) * ri for ri in rdir)
    maxdeg = sum(box)

    # Homogeneous-layer recurrence for y = 1/Q: writing Q = c0 + (higher order),
    # the degree-m part of y is y_m = -(1/c0) * sum_{|e|>=1} Q_e * y_{m-|e|}(shifted by e).
    layers = {0: {zero: 1 / c0}}
    Qterms = [(e, c, sum(e)) for e, c in Qd.items()]
    for m in range(1, maxdeg + 1):
        layer = {}
        for e, c, de in Qterms:
            src = layers.get(m - de)
            if not src:
                continue
            for ey, cy in src.items():
                enew = tuple(a + b for a, b in zip(ey, e))
                if any(a > b for a, b in zip(enew, box)):
                    continue
                layer[enew] = layer.get(enew, QQ(0)) - c * cy / c0
        if layer:
            layers[m] = layer
    ycoeff = {}
    for lay in layers.values():
        ycoeff.update(lay)

    # Multiply by P and extract the diagonal coefficients: the coefficient of
    # x^(j*rdir) in P*y only involves the finitely many terms of P.
    Pd = {tuple(e): QQ(c) for e, c in P.dict().items()}
    terms = []
    for j in range(N):
        target = tuple(j * ri for ri in rdir)
        s = QQ(0)
        for e, c in Pd.items():
            rem = tuple(a - b for a, b in zip(target, e))
            if any(a < 0 for a in rem):
                continue
            v = ycoeff.get(rem)
            if v is not None:
                s += c * v
        terms.append(s)
    return terms


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