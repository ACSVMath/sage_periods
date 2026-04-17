r"""
Methods for doing evaluation-interpolation, rational reconstruction, and the Chinese remainder theorem.
"""

from sage.arith.misc import CRT_list
from sage.rings.finite_rings.integer_mod_ring import Integers
from sage.matrix.constructor import Matrix
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.structure.sequence import Sequence
from sage.rings.integer_ring import ZZ
from sage.arith.misc import crt
from sage.arith.misc import gcd
from sage.arith.functions import lcm
from sage.matrix.constructor import matrix
from sage.structure.element import parent
from sage.misc.misc_c import prod
from sage.arith.misc import rational_reconstruction
from sage.modules.free_module_element import vector

from collections import Counter # For Stage 3

from sage.misc.verbose import verbose, set_verbose

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
      Raises ``ValueError`` if the Cor. 5.18 gcd test from von zur Gathen--Gerhard fails.

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

    # Return the leading coefficient, or 1 when input is zero.
    def lu(h):
        return h.leading_coefficient() if h != 0 else F.one()

    # Return the monic associate, or 0 when input is zero.
    def normal(h):
        return h / lu(h) if h != 0 else R.zero()

    # By remark below Cor 5.18 in vzG-G, to check that gcd(r_j,t_j) = 1 it suffices to
    # check that gcd(m,t_j) = 1. There is a solution if and only if this holds. 
    def admissible(den):
        if den != 0 and gcd(den,m) == 1:
            return True
        else:
            return False

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

    # Extended Euclidean Algorithm (EEA) initialization (with monic remainders)
    r0 = normal(m)
    t0 = R.zero()

    r1 = normal(g)
    t1 = R(lu(g)**(-1))   # So r1 = t1*g

    best = None
    if admissible(t1):
        best = (score(r1, t1), r1, t1)

    # We now run our EEA loop. Modifies things in-place, so
    # r2 = r_{i+1}, r1 = ri, r0 = r_{i-1}
    # t2 = t_{i+1}, t1 = ti, t0 = t_{i-1}.

    # For scoring, we run through the entire EEA keeping track
    # of the lowest-scored admissible pair that we encounter.
    while r1 != 0:
        q, rem = r0.quo_rem(r1)
        inv_rho = lu(rem)**(-1)
        r2 = inv_rho * rem
        t2 = inv_rho * (t0 - q * t1)

        if admissible(t2):
            cand = (score(r2, t2), r2, t2)
            if best is None or cand[0] < best[0]:
                best = cand

        r0, r1 = r1, r2
        t0, t1 = t1, t2

    if best is None:
        raise ValueError("No admissible Cauchy reconstruction found from these sample points.")

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
    num: int = 0
    points: set = field(default_factory=set)
    rand: dict = field(default_factory=dict)
    counter: int = 1
    nbticks: int = 1

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
        self.num=0
        self.points = set()
        self.rand = {}
        self.counter = 1
        self.nbticks=1
        return

    def reset(self):
        r"""
        Reset the reconstruction state to its empty initial configuration.

        This clears all cached samples, test values, current reconstructions,
        recorded points, and any candidate accumulated so far.

        INPUT:

        * ``None`` -- This method takes no user-supplied arguments.

        OUTPUT:

        * ``None`` -- This resets the ``ReconstructionData`` instance in place.
        """
        self.vals = {}
        self.recons = {}
        self.tests = {}
        self.candidate = None
        self.num = 0
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
            sage: # recon_add_rat only tries a full reconstruction after 1, 6, 15, ... samples.
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
    L_rho0, shape_rho0 = Sequence(val[0].list()), val[0].dimensions()
    L_B, shape_B = Sequence(val[1].list()), val[1].dimensions()

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

    # Do not try reconstruction at every new point: back off geometrically to keep the overall cost low
    H.nbticks -= 1
    if H.nbticks > 0:
        return
    H.counter += 4
    H.nbticks = H.counter

    # Get the points associated with these vals / tests
    points = Sequence(v[1] for v in all_tests)

    # Reconstruct the random linear form as a rational function of the parameter
    recon = cauchy_interp_rational_fn(K, zip([Kev(p) for p in points], [v[0] for v in all_tests]))

    # Accept the reconstruction only after obtaining the same candidate twice
    if key in H.recons and H.recons[key] == recon:
        if H.num > -1:
            all_vals = H.vals[key]
            cand = Sequence([])
            den = K(0)

            # First recover a common denominator from the same random linear
            # combination.  This is much cheaper than reconstructing denominators
            # for all coordinates separately.
            den = cauchy_interp_rational_fn(K, zip(points, [sum(all_vals[i][0][j] * rand[j] for j in range(len(L))) for i in range(len(points))]), criterion="total_degree").denominator()

            # Once the denominator is known, each numerator is just an ordinary interpolation problem.
            if den.is_constant():
                evden = [den for p in points]
            else:
                evden = [den(p) for p in points]
            
            '''
                TODO:
                Potential improvement: Cached batch reconstruction of the numerators? Or use one of the SOTA methods?
                And also cache our reconstruction of the denominator.
            '''
            for j in range(len(L)):
                cand.append(Rt.lagrange_polynomial(list(zip(points, Sequence(evden[i] * all_vals[i][0][j] for i in range(len(points)))))))

            # Now, cand is a sequence of same length as L, containing our reconstruction attempt
            if max((p.degree() for p in cand), default=ZZ(-1)) > QQ(2) * len(points) / 3:
                # Bad luck: the sampling is not good
                # This happens so rarely that we just throw away everything
                H.reset()
                return 

            # Otherwise, we have a good reconstruction candidate
            # Rebuild our pair entries separately, then assign candidate
            rho0_cand = Sequence(p / den for p in cand[:bkpt])
            B_cand = Sequence(p / den for p in cand[bkpt:])

            H.candidate = ( matrix(rho0_cand[0].parent(), *shape_rho0, rho0_cand), matrix(B_cand[0].parent(), *shape_B, B_cand) )
        else:
            H.num += 1

    H.recons[key] = recon
    return


##### Code for computing the dependency among the rho[i] #####

def compute_reductions_dependency(rho0,B):
    r"""
    Generates and computes a linear dependency among the $\rho_i$ as described in Section 7 of Lairez 2016.

    The matrix $\rho_0$ represents the reconstructed $\rho_0(t)$ with respect to a list of monomials $M$
    in $GF(p)(t)[x_0...x_n,u,v]$, all of which have the form $p(x)u^{n+1}$ (i.e., represent top level forms). 
    The matrix $B$ defines the map $m$ with respect to the basis ``M``.

    This is an internal function for sage_periods, and is not meant to be called by the user.

    INPUT:

    * ``rho0`` -- A Sage matrix over $GF(p)(t)$ representing the reconstructed class $\rho_0(t)$ in the basis $M$.
    * ``B`` -- A square Sage matrix over $GF(p)(t)$ representing the connection action on that basis.

    OUTPUT:

    * ``Lp`` -- A Python list of polynomials in $GF(p)[t]$ giving a nontrivial linear dependency among the iterates
      $\rho_0, \rho_1, \ldots$, with denominators cleared so that the result can be lifted across primes by CRT.

    EXAMPLES:

    A nontrivial example.

            sage: K = PolynomialRing(GF(7), "t").fraction_field()
            sage: t = K.gen()
            sage: rho0 = matrix(K, [[1]])
            sage: B = matrix(K, [[1/(t + 1)]])
            sage: # rho_1 = (t+1)^(-1) rho_0, so clearing denominators should give [1, t + 1]
            sage: compute_reductions_dependency(rho0, B)
            [1, t + 1]
            
    A trivial example.

            sage: K = PolynomialRing(GF(7), "t").fraction_field()
            sage: rho0 = matrix(K, [[1]])
            sage: B = matrix(K, [[0]])
            sage: # rho_1 = 0, so we should get [0,1] as output
            sage: compute_reductions_dependency(rho0, B)
            [0, 1]

    """
    K = B[0,0].parent() # K is F_p(t)
    t = K.gen() # This is the generator for GF(p)(t)

    m=0
    # All our rho[i]'s should be represented by vectors with respect to M
    # Since rho0 is a matrix, this is easy
    rho = [vector(rho0.list())]
    while True:
        verbose("        m: "+str(m),level=1)
        verbose("rho: ",level=1)
        verbose(rho,level=1)
        # Note that our degree is guaranteed to stay bounded (i.e. r need not increase)
        # because the M we found in Stage 1 is canonical, as are rho0 and B.
        rhomat = Matrix(rho).transpose()
        
        if rhomat.rank() == m + 1:
            rho.append(vector([ fn.derivative(t)for fn in rho[m]]) + B*rho[m])
        else:
            # Need to solve the system over K to find coefficients for 
            # which the linear comb of the lower rho_k = rho_m
            try:
                y = rhomat[:,:m].solve_right(rho[m])
            except:
                print("            OH NO! We couldn't find a dependency!")
            
            # At this point, we have our coefficient vector. We could return, but we should clear denominators first to make CRT possible.

            # Extract denominators of y coefficients
            denoms = [ai.denominator() for ai in y]
            Rp = K.ring() # Rp = GF(p)[t] and each entry in denoms is already in Rp
            # Take LCM of denoms
            LCM = lcm(denoms) if denoms else Rp(1) # If our denoms list is empty then our LCM is 1
            out = [LCM*y[i] for i in range(m)] + [LCM] # This could have entries in K or Rp
            try:
                out = [Rp(c) for c in out]
            except Exception:
                raise Exception("BAD_PRIME")
            return out
        m += 1

def lift_operator_across_primes(t, Lp_coeffs_dict, bad_primes):
    r"""
    Reconstruct operator coefficients in $\mathbb{Q}[t]$ from reductions in $GF(p)[t]$
    using CRT and rational reconstruction.

    This is an internal function for sage_periods, and is not meant to be called by the user.

    INPUT:

    * ``t`` -- The generator that we should also rebuild our operators in
      (``t`` should match across the entries of ``Lp_coeffs_dict``).
    * ``Lp_coeffs_dict`` -- A Python dictionary whose keys are primes and whose values are the Stage 2 outputs,
      consisting of lists of polynomials constituting an annihilating operator modulo primes $p$.
    * ``bad_primes`` -- A Python list of primes that we should avoid considering in
      our reconstructions. We both read from and write to this list in this function.

    OUTPUT:

    * ``L`` -- A Python list of coefficients in $\mathbb{Q}[t]$ obtained by lifting the modular operators.

    EXAMPLES:

    A non-trivial example.

            sage: RQQ = PolynomialRing(QQ, "t")
            sage: t = RQQ.gen()
            sage: R5 = PolynomialRing(GF(5), "t")
            sage: t5 = R5.gen()
            sage: R7 = PolynomialRing(GF(7), "t")
            sage: t7 = R7.gen()
            sage: R11 = PolynomialRing(GF(11), "t")
            sage: t11 = R11.gen()
            sage: Lp_coeffs_dict = {
            sage:     5: [3*t5 + 2, t5**2 + 1],
            sage:     7: [4*t7 + 5, t7**2 + 1],
            sage:     11: [6*t11 + 4, t11**2 + 1],
            sage: }
            sage: bad_primes = []
            sage: lift_operator_across_primes(t, Lp_coeffs_dict, bad_primes)
            [1/2*t + 1/3, t^2 + 1]
            sage: bad_primes
            []
            
    A trivial example.

            sage: RQQ = PolynomialRing(QQ, "t")
            sage: t = RQQ.gen()
            sage: R5 = PolynomialRing(GF(5), "t")
            sage: R7 = PolynomialRing(GF(7), "t")
            sage: Lp_coeffs_dict = {
            sage:     5: [R5(1)],
            sage:     7: [R7(1)],
            sage: }
            sage: bad_primes = []
            sage: lift_operator_across_primes(t, Lp_coeffs_dict, bad_primes)
            [1]
            sage: bad_primes
            []
            
    An example with a bad prime.

            sage: RQQ = PolynomialRing(QQ, "t")
            sage: t = RQQ.gen()
            sage: R5 = PolynomialRing(GF(5), "t")
            sage: t5 = R5.gen()
            sage: R7 = PolynomialRing(GF(7), "t")
            sage: t7 = R7.gen()
            sage: R11 = PolynomialRing(GF(11), "t")
            sage: t11 = R11.gen()
            sage: R13 = PolynomialRing(GF(13), "t")
            sage: t13 = R13.gen()
            sage: Lp_coeffs_dict = {
            sage:     5: [3*t5 + 2, t5**2 + 1],
            sage:     7: [4*t7 + 5, t7**2 + 1],
            sage:     11: [6*t11 + 4, t11**2 + 1],
            sage:     13: [7*t13 + 9, t13 + 1],
            sage: }
            sage: bad_primes = []
            sage: lift_operator_across_primes(t, Lp_coeffs_dict, bad_primes)
            [1/2*t + 1/3, t^2 + 1]
            sage: bad_primes
            [13]
    """

    # Avoid duplicate bad primes
    def _append_bad(p):
        if p not in bad_primes:
            bad_primes.append(p)

    # Scale coefficients so that leading coefficient (in t) of final entry is 1 in GF(p)
    def _normalize(coeffs):
        if not coeffs:
            return None
        lead_poly = coeffs[-1]
        if lead_poly == 0:
            return None
        lc = lead_poly.lc()
        if lc == 0:
            return None
        scale = lc**(-1)
        return [scale * c for c in coeffs]

    # Handle trivial case
    if not Lp_coeffs_dict:
        verbose("        Lp_coeffs_dict is empty. Need more primes...",level=1)
        raise Exception("Lp_coeffs_dict is empty. Need more primes...")

    #### Pre-CRT pruning: normalization + consistency checks ####
    # We use heuristics to make sure that our operators are all the same "shape," i.e. non-degenerate.
    changed = True
    while changed and Lp_coeffs_dict:
        changed = False
        to_drop = set()

        # Normalize each L_p
        for p in list(Lp_coeffs_dict.keys()):
            try:
                scaled = _normalize(Lp_coeffs_dict[p])
            except Exception:
                to_drop.add(p)
                verbose(f"Dropped prime {p} because of failed normalization. This means that either our coeff list was empty or had zero leading term.",level=1)
            else:
                Lp_coeffs_dict[p] = scaled

        if to_drop:
            changed = True
            for p in to_drop:
                Lp_coeffs_dict.pop(p, None)
                _append_bad(p)
            continue

        # Enforce common length in coefficient list
        lengths = {p: len(Lp_coeffs_dict[p]) for p in Lp_coeffs_dict}
        mode_len = Counter(lengths.values()).most_common(1)[0][0]
        for p, L in lengths.items():
            if L != mode_len:
                to_drop.add(p)
                verbose(f"Dropped prime {p} because the coefficient list was not of the correct length.",level=1)

        if to_drop:
            changed = True
            for p in to_drop:
                Lp_coeffs_dict.pop(p, None)
                _append_bad(p)
            continue

        # Enforce common degrees per coefficient polynomial
        primes = sorted(Lp_coeffs_dict.keys())
        mplus1 = mode_len
        for i in range(mplus1):
            degs = {p: Lp_coeffs_dict[p][i].degree() for p in primes}
            expected_deg = Counter(degs.values()).most_common(1)[0][0]
            for p, d in degs.items():
                if d != expected_deg:
                    to_drop.add(p)
                    verbose(f"Dropped prime {p} because one of its coefficients has degree degeneration, so a lift won't work.",level=1)
        if to_drop:
            changed = True
            for p in to_drop:
                Lp_coeffs_dict.pop(p, None)
                _append_bad(p)

    if not Lp_coeffs_dict:
        verbose("        No good primes left after pruning. Need more primes...",level=1)
        raise Exception("No good primes left after pruning. Need more primes...")

    verbose("        Pruned bad primes. Performing CRT + reconstruction...",level=1)
    verbose("        FINAL Lp_coeffs_dict after pruning, canonical entries: ",level=1)
    verbose(Lp_coeffs_dict,level=1)
    primes = sorted(Lp_coeffs_dict.keys())
    mplus1 = len(next(iter(Lp_coeffs_dict.values()))) # These are all common length now
    M = ZZ(prod(primes))

    # Build QQ[t]
    RQQ = PolynomialRing(QQ, str(t))
    t = RQQ.gen()

    #### CRT + rational reconstruction ####
    '''
        TODO: 
        - Use Sage's CRT_vectors() call, which takes in a list of integer vectors
          and does vector-wise CRT on that. Faster than doing it per-entry.
        - Use the MultiModularBasis class, which does incremental CRT
          but only over the integers and only a single entry at a time.
          Faster than "naive" CRT adding points one at a time.
    '''
    lifted = []
    for i in range(mplus1): 
    # This layer represents the coefficients a_i in our final operator
        deg_i = Lp_coeffs_dict[primes[0]][i].degree() # This is common across mods
        if deg_i < 0:
            lifted.append(RQQ.zero())
            continue

        coeffs_i = []
        for k in range(deg_i + 1): 
            # This layer represents the terms in a_i
            residues = []
            for p in primes: # Reconstruct that term across primes
                ai_p = Lp_coeffs_dict[p][i]
                ck = ai_p[k]
                residues.append(ZZ(ck))
            x = ZZ(CRT_list(residues, primes))  # 0 <= x < M

            # Use rational reconstructiont to find a/b congruent to x (mod M)
            try:
                q = rational_reconstruction(x, M)  # May raise ValueError
                coeffs_i.append(QQ(q))
            except (ValueError, ArithmeticError):
                verbose("        Rational reconstruction failed for some coefficient. Need more primes...",level=1)
                raise Exception(
                    "Rational reconstruction failed for some coefficient. Need more primes..."
                )

        lifted.append(RQQ(coeffs_i))
    verbose("        Successfully lifted! :) Returning operator.",level=1)
        
    return lifted