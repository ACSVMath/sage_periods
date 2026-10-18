r"""
Methods for doing evaluation-interpolation, rational reconstruction, and the Chinese remainder theorem.
"""

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.rings.integer_ring import ZZ
from sage.arith.misc import gcd
from sage.matrix.constructor import matrix
from sage.structure.element import parent
from sage.misc.misc_c import prod

from .errors import BadPrimeError, NeedMorePrimesError

# For slots in ReconstructionData class
from dataclasses import dataclass, field


##### Code for Cauchy interpolation #####

def cauchy_interp_rational_fn(K, pairs, criterion="total_degree"):
    r"""
    Performs Cauchy interpolation across all k from 0 to len(pairs) inclusive and chooses
    the "best" rational output according to ``criterion``.

    The possible criterion are

    - ``total_degree``: minimize ``deg(num)+deg(den)``
    - ``max_degree``:   minimize ``max(deg(num),deg(den))``
    - ``den_degree``:   minimize ``deg(den)`` then ``deg(num)``

    The Cauchy interpolation occurs over $\mathbb{F}_p$, following von zur Gathen--Gerhard Section 5.8, Corollary 5.18.

    This is an internal function for sage_periods, and is not meant to be called by the user.

    INPUT:

    * ``K`` -- A Sage fraction field, typically $GF(p)(t)$, in which the reconstructed function $h(t)$ lies.
    * ``pairs`` -- A Python list ``[(u_i, v_i), ...]`` with each $v_i$ in $GF(p)$ and each $u_i$ coercible to $GF(p)$
      (for example $u_i \in \mathbb Q$ with denominator not divisible by $p$).
    * ``criterion`` -- (Optional) A string specifying how to choose the "best" admissible reconstruction among the EEA candidates.
      Supported values are ``total_degree``, ``max_degree``, and ``den_degree``.

    OUTPUT:

    * ``h`` -- An element of $GF(p)(t)$ such that $h(u_i) = v_i$ for all sample pairs,
      with $\deg(\mathrm{numer}) < k$ and $\deg(\mathrm{denom}) \leq \mathrm{len}(\mathrm{pairs})-k$ for some admissible $k$ if the reconstruction is solvable.
      The gcd test from Cor. 5.18 excludes denominators vanishing at sample points.
      Polynomial interpolation (denominator 1) always supplies an admissible candidate.

    ASSUMPTIONS: (not checked or handled gracefully)

    * None of the points $u_i$ have denominators divisible by $p$.
    * The $u_i$ are distinct.
    * That $0 \leq k < \mathrm{len}(\mathrm{pairs})$.

    EXAMPLES:

    A nontrivial example.

            sage: K = PolynomialRing(GF(97),"t").fraction_field()
            sage: t = K.gen()
            sage: Fp = K.base_ring()
            sage: f = K.random_element()
            sage: f
            (70*t^2 + 38*t + 37)/(t^2 + 60*t + 8)
            sage: u = [ Fp(i) for i in range(20)] 
            sage: v = [f(u[i]) for i in range(20)]
            sage: pairs = list(zip(u,v))
            sage: cauchy_interp_rational_fn(K,pairs)
            (70*t^2 + 38*t + 37)/(t^2 + 60*t + 8)        

    A trivial example.

            sage: v = [Fp(0) for i in range(20)]
            sage: pairs = list(zip(u,v))
            sage: cauchy_interp_rational_fn(K,pairs)
            0
            
    Changing criteria can be useful, depending on the application.

            sage: K = PolynomialRing(GF(97),"t").fraction_field()
            sage: t = K.gen()
            sage: Fp = K.base_ring()
            sage: f = K.random_element()
            sage: f
            (77*t^2 + 59*t + 12)/(t^2 + 36*t + 1)
            sage: u = [ Fp(i) for i in range(20)] # This will probably be larger than the degree of f.
            sage: v = [f(u[i]) for i in range(20)]
            sage: pairs = list(zip(u,v))
            sage: cauchy_interp_rational_fn(K,pairs,criterion='total_degree')
            (77*t^2 + 59*t + 12)/(t^2 + 36*t + 1)
            sage: cauchy_interp_rational_fn(K,pairs,criterion='den_degree')
            72*t^19 + 27*t^18 + 58*t^17 + 39*t^16 + 41*t^15 + 52*t^14 + 81*t^13 + 92*t^12 + 72*t^11 + 30*t^10 + 70*t^9 + 29*t^8 + 21*t^7 + 51*t^6 + 22*t^5 + 79*t^4 + 40*t^3 + 41*t^2 + 50*t + 12

    """

    F = K.base_ring() # GF(p)
    R = K.ring() # GF(p)[t]
    t = R.gen()

    pts = [(F(u), F(v)) for (u, v) in pairs]
    g = R.lagrange_polynomial(pts)
    m = prod((t - ui) for (ui, _) in pts)

    # Evaluate a given (numerator, denominator) pair based on input criterion.
    def score(num, den):
        dn = num.degree()
        dd = den.degree()
        if criterion == "total_degree":
            return (dn + dd, dd, dn)
        if criterion == "max_degree":
            return (max(dn, dd), dd, dn)
        if criterion == "den_degree":
            return (dd, dn)
        raise ValueError("Unknown criterion")

    # Ordinary EEA suffices: scaling a remainder and its Bezout coefficient
    # together does not change their quotient, degrees, or admissibility.
    r0, r1 = m, g
    t0, t1 = R.zero(), R.one()
    best = (score(r1, t1), r1, t1)  # Polynomial interpolation is admissible.

    while r1 != 0:
        q, r2 = r0.quo_rem(r1)
        t2 = t0 - q * t1
        candidate_score = score(r2, t2)
        # A gcd is only needed for candidates that can improve the score.
        if candidate_score < best[0] and t2 != 0 and gcd(t2, m) == 1:
            best = (candidate_score, r2, t2)
        r0, r1 = r1, r2
        t0, t1 = t1, t2

    _, num_best, den_best = best

    # Normalize denominator to make it monic
    lc = den_best.leading_coefficient()
    if lc != 1:
        inv_lc = lc**(-1)
        num_best *= inv_lc
        den_best *= inv_lc

    return K(num_best) / K(den_best)

'''
    TODO:
    Potential improvement to coordinatewise Cauchy Interpolation by caching, or using
    asymptotically faster methods for vector rational reconstruction.
'''


@dataclass(slots=True)
class ReconstructionData:
    r"""
    Container storing the mutable state used by evaluation/interpolation and CRT reconstruction.

    This record keeps exact samples, compressed test values, current rational reconstructions,
    the set of already used points, and a candidate once reconstruction has stabilized.
    """
    vals: dict = field(default_factory=dict)
    recons: dict = field(default_factory=dict)
    tests: dict = field(default_factory=dict)
    candidate: object = None
    points: set = field(default_factory=set)
    rand: dict = field(default_factory=dict)
    counter: int = 1
    nbticks: int = 1
    first_tick: object = None  # Optional sample-count hint from an earlier prime.

    def __init__(self):
        r"""Create an empty reconstruction state record.

        The record stores the data accumulated during CRT lifting or
        evaluation/interpolation: exact samples, compressed test values, current
        reconstructions, the set of already used points, and a candidate once
        reconstruction has stabilized.

        INPUT:

        * ``None`` -- This initializer takes no user-supplied arguments.

        OUTPUT:

        * ``None`` -- This initializes the ``ReconstructionData`` instance in place.
        """
        self.vals = {}
        self.recons = {}
        self.tests = {}
        self.candidate = None
        self.points = set()
        self.rand = {}
        self.counter = 1
        self.nbticks=1
        self.first_tick = None
        return

    def reset(self):
        r"""
        Reset the reconstruction state to its empty initial configuration.

        This clears all cached samples, test values, current reconstructions,
        recorded points, and any candidate accumulated so far. The scheduling
        hint ``first_tick`` is preserved; it is not reconstructed data.

        INPUT:

        * ``None`` -- This method takes no user-supplied arguments.

        OUTPUT:

        * ``None`` -- This resets the ``ReconstructionData`` instance in place.
        """
        self.vals = {}
        self.recons = {}
        self.tests = {}
        self.candidate = None
        self.points = set()
        self.rand = {}
        self.counter = 1
        self.nbticks = 1    

def recon_add_rat(H,val,point,xkey):
    r"""
    Reconstruction manager for evaluation-interpolation in our Gauss-Manin stage.

    Here, ``val`` is a pair of the form $(\rho_0, B)$ where $\rho_0$ and $B$ are evaluated
    objects obtained from ``gauss_manin_helper``. They're both matrices, but $\rho_0$ has dimensions 
    $|M| \times 1$ matrix while $B$ has dimensions $|M| \times |M|$.

    This is an internal function for sage_periods, and is not meant to be called by the user.

    INPUT:

    * ``H`` -- A ``ReconstructionData`` instance storing the state of the current evaluation/interpolation reconstruction.
    * ``val`` -- A pair ``(rho0, B)`` of evaluated Sage matrices, expressed in the same sampled basis.
    * ``point`` -- The evaluation point of the parameter at which ``val`` was computed.
    * ``xkey`` -- A hashable key encoding the sampled basis pattern; only samples with the same key are combined.

    OUTPUT:

    * ``None`` -- This function updates ``H`` in place by adding the new sample and, when enough consistent data is available,
      storing a reconstructed candidate in ``H.candidate``.

    EXAMPLES:

    Nontrivial example:

            sage: K = PolynomialRing(GF(101), "t").fraction_field()
            sage: t = K.gen()
            sage: F = K.base_ring()
            sage: H = ReconstructionData()
            sage: rho0_fun = matrix(K, [[t + 1], [2*t + 3]])
            sage: B_fun = matrix(K, [[1, t], [t**2, t + 2]])
            sage: def eval_matrix(M, u):
            sage:     return matrix(F, M.nrows(), M.ncols(), [entry(u) for entry in M.list()])
            sage: 
            sage: H.candidate
            None
            sage: # recon_add_rat only tries a full reconstruction after 1, 3, 6, 10, 15, ... samples.
            sage: # Using 15 points lets the stabilized reconstruction appear twice.
            sage: for u in range(15):
            sage:     u_F = F(u)
            sage:     rho0_eval = eval_matrix(rho0_fun, u_F)
            sage:     B_eval = eval_matrix(B_fun, u_F)
            sage:     recon_add_rat(H, (rho0_eval, B_eval), K(u), "fixed basis")
            sage: 
            sage: H.candidate
            ([  t + 1]
            [2*t + 3], [    1     t]
            [  t^2 t + 2])
            
    A trivial example.

            sage: K = PolynomialRing(GF(5), "t").fraction_field()
            sage: F = K.base_ring()
            sage: H = ReconstructionData()
            sage: rho0_eval = matrix(F, 0, 1, [])
            sage: B_eval = matrix(F, 0, 0, [])
            sage: H.candidate
            None
            sage: recon_add_rat(H, (rho0_eval, B_eval), K(0), "empty basis")
            sage: H.candidate
            ([], [])

    """

    # Flatten each matrix in val
    L_rho0, shape_rho0 = val[0].list(), val[0].dimensions()
    L_B, shape_B = val[1].list(), val[1].dimensions()

    # Handle corner cases
    if len(L_rho0) == 0 and len(L_B) == 0:
        H.points.add(point)
        H.candidate = (val[0], val[1])
        return
    if len(L_rho0) == 0 or len(L_B) == 0:
        raise ValueError("Inconsistent empty stage-1 data")

    # Keep track of "breakpoint" between these two objects in conjoined vector,
    # i.e. the length of first vector, which is the index of first entry of second vector.
    bkpt = len(L_rho0)

    # Join the two vectors into one; we'll run our reconstruction machinery on this
    # single vector, then separate and rebuild the pair at the end using bkpt.
    assert parent(L_rho0[0]) == parent(L_B[0]), "rho0 and B evaluations must have same universe!"
    L = L_rho0 + L_B
    Kev = L[0].parent()
    K = point.parent()
    Rt = K.ring()
    
    key = xkey # M' exponent basis
    H.points.add(point)

    # Following code of Lairez, we use a fixed random linear form on the flattened coordinates.  This gives a
    # cheap scalar test for whether the current sample set is compatible with a stable rational reconstruction.
    if key in H.rand:
        rand = H.rand[key]
    else:
        rand = [Kev.random_element() for _ in L]
        H.rand[key] = rand

    # Compress the whole object to one scalar-valued sample used only as a
    # consistency test before the more expensive full coefficient reconstruction
    test = sum((rand[i] * L[i] for i in range(len(L))), Kev.zero())

    # Add this test and our value to "buckets"
    if key in H.vals:
        H.vals[key].append((L,point))
    else:
        H.vals[key] = [(L,point)]
    if key in H.tests:
        H.tests[key].append((test,point))
    else:
        H.tests[key] = [(test,point)]

    all_tests = H.tests[key]

    # A hint schedules attempts, not acceptance: stability is still checked below.
    if H.first_tick is not None:
        if len(all_tests) < H.first_tick:
            return
    else:
        # Without a hint, try after 1, 3, 6, 10, ... samples.
        H.nbticks -= 1
        if H.nbticks > 0:
            return
        H.counter += 1
        H.nbticks = H.counter

    # Interpolation data live in GF(p), not in the larger fraction field.
    # In particular, den(u) should be a field element, not a constant fraction.
    points = [Kev(v[1]) for v in all_tests]

    # Reconstruct the random linear form as a rational function of the parameter
    recon = cauchy_interp_rational_fn(K, zip(points, [v[0] for v in all_tests]))

    # Accept the reconstruction only after obtaining the same candidate twice
    if key in H.recons and H.recons[key] == recon:
        all_vals = H.vals[key]
        cand = []

        # Reuse the denominator of the same random linear combination.
        den = recon.denominator()

        # Once the denominator is known, each numerator is just an ordinary interpolation problem.
        if den.is_constant():
            evden = [Kev(den[0])] * len(points)
        else:
            evden = [den(p) for p in points]
        
        '''
            TODO:
            Potential improvement: Cached batch reconstruction of the numerators? Or use one of the SOTA methods?
        '''
        for j in range(len(L)):
            cand.append(Rt.lagrange_polynomial([(u, evden[i] * all_vals[i][0][j]) for i, u in enumerate(points)]))

        # Now, cand is a sequence of same length as L, containing our reconstruction attempt
        if max((p.degree() for p in cand), default=ZZ(-1)) > QQ(2) * len(points) / 3:
            # Bad luck: the sampling is not good
            # This happens so rarely that we just throw away everything.
            # Also drop the hint, so that we fall back to the triangular
            # schedule instead of retrying at the same too-small count forever.
            H.reset()
            H.first_tick = None
            return 

        # Otherwise, we have a good reconstruction candidate
        # Rebuild our pair entries separately, then assign candidate
        rho0_cand = [p / den for p in cand[:bkpt]]
        B_cand = [p / den for p in cand[bkpt:]]

        H.candidate = ( matrix(rho0_cand[0].parent(), *shape_rho0, rho0_cand), matrix(B_cand[0].parent(), *shape_B, B_cand) )

    H.recons[key] = recon
    return


##### Incremental CRT lifting across primes (Stage 3) #####

class CRTLifter:
    r"""
    Incremental Chinese-remainder accumulator for Stage 3.

    Residue vectors over $GF(p)$ (one list of integers per prime) are folded in as they arrive:
    each fold costs a single CRT step per entry, and each lift attempt is one rational
    reconstruction pass over the accumulated integers. Vectors are keyed by a hashable
    ``shape`` (e.g. a degree pattern), and only primes with identical shapes are combined.
    A prime with a degenerate shape accumulates in a minority class: it never corrupts
    the majority, and it is kept (not discarded), in case its class later becomes the majority.

    This is an internal class for sage_periods, and is not meant to be used by the user.

    EXAMPLES:

    A wrong-shape prime folded first does not affect the majority lift.

            sage: from sage_periods.reconstruction import CRTLifter
            sage: L = CRTLifter()
            sage: L.fold(5, [1, 2], ("bad",))
            sage: for p in [7, 11, 13]:
            ....:     L.fold(p, [ZZ(1/2 % p), ZZ(-1/3 % p), 1], ("good",))
            sage: L.lift_majority()
            (('good',), [1/2, -1/3, 1], frozenset({7, 11, 13}))
            sage: L.class_of(5), L.class_of(7), L.class_of(3)
            (('bad',), ('good',), None)
    """

    def __init__(self):
        r"""Create an empty accumulator."""
        self.classes = {}   # shape -> {"modulus", "residues", "primes"}
        self.shape_of = {}  # prime -> shape of the class it was folded into

    def fold(self, p, residues, shape):
        r"""
        Fold the residues ``residues`` (integers modulo ``p``) into the class ``shape``.
        Folding a prime twice has no effect.
        """
        p = ZZ(p)
        if p in self.shape_of:
            return
        self.shape_of[p] = shape
        c = self.classes.get(shape)
        if c is None:
            self.classes[shape] = {
                "modulus": p,
                "residues": [ZZ(a) % p for a in residues],
                "primes": {p},
            }
            return
        M = c["modulus"]
        Minv = M.inverse_mod(p)
        # x' = x + M * ((a - x) * M^{-1} mod p) is the unique lift with
        # x' = x mod M and x' = a mod p, 0 <= x' < M p.
        c["residues"] = [x + M * (((ZZ(a) - x) * Minv) % p)
                         for x, a in zip(c["residues"], residues)]
        c["modulus"] = M * p
        c["primes"].add(p)

    def _majority(self):
        r"""Return ``(shape, class)`` for the class with the most primes (first-inserted on ties)."""
        return max(self.classes.items(), key=lambda kv: len(kv[1]["primes"]))

    def lift_majority(self):
        r"""
        Rationally reconstruct the class with the most primes.

        OUTPUT:

        * ``(shape, values, primes)`` with ``values`` a list of rationals and ``primes``
          the frozenset of primes in that class. Raises ``NeedMorePrimesError`` if no prime
          has been folded yet or if rational reconstruction fails for some entry.
        """
        if not self.classes:
            raise NeedMorePrimesError("No accumulated primes. Need more primes...")
        shape, c = self._majority()
        M = c["modulus"]
        out = []
        for x in c["residues"]:
            try:
                out.append(x.rational_reconstruction(M))
            except (ValueError, ArithmeticError):
                raise NeedMorePrimesError(
                    "Rational reconstruction failed for some coefficient. Need more primes..."
                )
        return shape, out, frozenset(c["primes"])

    def class_of(self, p):
        r"""Return the shape of the class that ``p`` was folded into, or ``None``."""
        return self.shape_of.get(ZZ(p))

    def majority_primes(self):
        r"""Return the frozenset of primes in the majority class (empty if nothing was folded)."""
        if not self.classes:
            return frozenset()
        return frozenset(self._majority()[1]["primes"])


##### Operator codec: residue vectors for CRTLifter #####

def operator_residues(coeffs):
    r"""
    Canonical residue vector of a modular operator, given as the Stage 2 coefficient list in $GF(p)[t]$.

    The list is scaled so that the leading polynomial is monic in $t$, then flattened densely by degree.
    Returns ``(shape, ints)`` with ``shape = (length, degree tuple)``. Raises ``BadPrimeError``
    if the list is empty or its leading polynomial is zero.

    This is an internal function for sage_periods, and is not meant to be called by the user.

    EXAMPLES:

            sage: from sage_periods.reconstruction import operator_residues
            sage: R7 = PolynomialRing(GF(7), "t"); t7 = R7.gen()
            sage: operator_residues([3*t7 + 1, 2*t7**2 + 4])
            ((2, (1, 2)), [4, 5, 2, 0, 1])
    """
    if not coeffs or coeffs[-1] == 0:
        prime = coeffs[-1].parent().characteristic() if coeffs else None
        raise BadPrimeError("Degenerate modular operator (empty or zero leading coefficient).", prime, "lifting")
    lc = coeffs[-1].leading_coefficient()
    scaled = [c / lc for c in coeffs]
    shape = (len(scaled), tuple(c.degree() for c in scaled))
    ints = []
    for c in scaled:
        ints.extend(ZZ(c[k]) for k in range(c.degree() + 1))
    return shape, ints


def operator_from_lift(t, shape, values):
    r"""
    Rebuild the $\mathbb{Q}[t]$ coefficient list from ``operator_residues`` data
    (``shape``, and the lifted rational ``values``).

    EXAMPLES:

            sage: from sage_periods.reconstruction import operator_from_lift
            sage: operator_from_lift("t", (3, (1, -1, 0)), [1/2, 1/3, 1])
            [1/3*t + 1/2, 0, 1]
    """
    RQQ = PolynomialRing(QQ, str(t))
    _, degs = shape
    out = []
    pos = 0
    for d in degs:
        if d < 0:
            out.append(RQQ.zero())
        else:
            out.append(RQQ(values[pos:pos + d + 1]))
            pos += d + 1
    return out


def operator_residue_count(shape):
    r"""Number of residues described by an ``operator_residues`` shape."""
    _, degs = shape
    return sum(d + 1 for d in degs if d >= 0)