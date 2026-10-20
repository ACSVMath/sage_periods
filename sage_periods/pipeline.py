r"""Main stages in our modular pipeline for the reduction algorithm"""

from sage.all import *
from .reconstruction import ReconstructionData, recon_add_rat
from .rham_koszul import RhamKoszulData, gauss_manin_helper
from .errors import *

from sage.misc.verbose import verbose, set_verbose

def prepare_coefficient_data(a, f):
    """Immutable characteristic-zero tables for f, fdelta, and a."""
    def table(coefficients):
        return tuple((e, c.numerator(), c.denominator())
                     for e, c in coefficients if c)
    fdict = f.dict()
    return (table(fdict.items()),
            table((e, c.derivative()) for e, c in fdict.items()),
            table(a.dict().items()))


# Stage 1
def compute_gauss_manin_connection(a,f,r,p,profile=None,degree_hint=None,staircase_keys=None,coefficient_data=None):
    r"""
    Compute the Gauss-Manin matrix $B$ for the map defined by closing $\rho_0$ under the map $\rho \mapsto [f^\delta \rho]_r$,
    as described in Section 7.2 of Lairez 2016. Also returns the basis for the space and the projection of $a$ into this space.

    This function is typically called as a subroutine of ``compute_period_annihilator``.

    INPUT:

    * ``a`` -- A polynomial representing the numerator of the function under consideration.
    * ``f`` -- A square-free polynomial representing the square-free denominator of the function under consideration.
    * ``r`` -- Order of the reduction $[\cdot]_r$ which we compute.
    * ``p`` -- The characteristic of the field which we evaluate ``a`` and ``f`` into to do our reductions.
    * ``profile`` -- (Optional) A ``ProfileData`` object shared across primes and points, modified in place by this function.
      implementing the "profile" variant of the MAGMA package for ``r = 2``.
      It recors which syzygy staircase lifts contribute to our reduction,
      so later points only enumerate those. Tremendous time-saver.
    
    * ``degree_hint`` -- (Optional) A sample count before the first reconstruction
      attempt; subsequent attempts occur at each new point. This is a scheduling
      hint only, and does not bypass the reconstruction stability check.

    * ``staircase_keys`` -- (Optional) A set recording staircase cache keys used by a probe.

    OUTPUT: 
    
    * $M$, $\rho_0(t)$ and $B(t)$ as described in Section 7.2 of Lairez 2016.

    """
    A = a.parent()
    QQ_t = A.base_ring().ring()
    bad_points = set()

    # Profiling applies only at r = 2: start recording if no profile is active yet.
    if profile is not None and r == 2 and profile.state == "disabled":
        profile.state = "recording"

    # Coefficient tables retain the parameter-polynomial layer GF(p)[t],
    # while evaluated forms belong to AF = GF(p)[x_0, ..., x_n].
    F = GF(p)
    Rp = PolynomialRing(F, QQ_t.variable_name())
    K = Rp.fraction_field()
    AF = A.change_ring(F) # This should equal RK.A across all RK

    # Reconstruction object for B and rho0
    H = ReconstructionData()
    H.first_tick = degree_hint

    if coefficient_data is None:
        coefficient_data = prepare_coefficient_data(a, f)

    def _build_coefficient_table(source):
        """Reduce the fixed coefficient data modulo this prime."""
        table = [(e, Rp(num), Rp(den)) for e, num, den in source]
        if any(den.is_zero() for _, _, den in table):
            raise BadPrimeError(
                "A coefficient denominator vanishes identically modulo p.",
                p,
                "gauss_manin",
            )
        return table

    def _evaluate_coefficient_table(table, pt):
        r"""
        Evaluate a coefficient table (exponent, numerator, denominator) at the
        point ``pt``, raising ``ZeroDivisionError`` when a coefficient
        denominator vanishes there.
        """
        out = {}
        for e, num, den in table:
            d = den(pt)
            if not d:
                raise ZeroDivisionError("coefficient denominator vanishes at this point")
            cv = num(pt) / d
            if cv:
                out[e] = cv
        return AF(out)

    try:
        ftab, fdeltatab, atab = map(_build_coefficient_table, coefficient_data)
    except BadPrimeError:
        raise
    except (ArithmeticError, TypeError, ValueError) as exc:
        raise BadPrimeError(
            "A coefficient of a, f, or f^delta cannot be reduced modulo p.",
            p,
            "gauss_manin",
        ) from exc

    # Counter for seeing if we must increase r, or if we just had bad evaluation point
    rtoosmall = 0

    # Main loop
    while True:
        u = ZZ.random_element(1, p)   # random evaluation point
        uF = F(u)
        u_K = K(uF)                   # for interpolation / storage mod p
        if u_K in H.points or uF in bad_points:
            continue
        verbose("        u: "+str(u),level=1)

        # Try evaluating f, a, f^delta into our ring. If this fails, pick a new point.
        try:
            feval = _evaluate_coefficient_table(ftab, uF)
            fdeltaeval = _evaluate_coefficient_table(fdeltatab, uF)
            aeval = _evaluate_coefficient_table(atab, uF)
        except ZeroDivisionError:
            verbose("Evaluations of f, a, or f^delta failed. Pick a different point.",level=1)
            bad_points.add(uF)
            continue
            
        # Build our RhamKoszulData object with evaluated f
        U = RhamKoszulData(feval,r=r,profile=profile,staircase_keys=staircase_keys)
        
        # Run gauss_manin_helper on <fteval,[aeval]>, with our prime and r
        try:
            ret = gauss_manin_helper(U, fdeltaeval, [aeval], None)
        except ReductionOrderTooSmallError:
            rtoosmall += 1
            bad_points.add(uF)
            if rtoosmall >= 3:
                raise ReductionOrderTooSmallError(
                    "Exceeded filtration bound at least 3 times; increase r.",
                    r,
                )
            continue

        # The point succeeded: offer its profile as a candidate. The profile
        # is armed once two consecutive candidates agree.
        if U.record_profile:
            profile.add_profile_candidate(U.profile_candidate(), p)

        # Extract results and add to interpolation routine. Proj is the |M| x 1 matrix expressing rho_0' in terms of M
        basis_key = ret.ebasis #gauss_manin_helper should return a tuple of tuples of tuples
        rho0_eval = ret.proj
        B_prime = ret.gm
        basis = ret.basis

        # Feed to interpolation manager
        recon_add_rat(H,(rho0_eval,B_prime),u_K,basis_key) # Both rho0 and B are matrices

        # Stability is handled by recon_add_rat, we just need to check if candidate has been assigned
        if H.candidate is not None:
            return H.candidate[0], basis, H.candidate[1] # rho0, M, B
        # Else, continue.

def _extend_echelon(basis, row):
    """Add a finite-field row if independent; inputs are left unchanged.

    Rows stay in discovery order. Each is zero at all earlier pivots, so
    eliminating in that order never reintroduces an already cleared pivot.
    """
    for pivot, previous in basis:
        if row[pivot]:
            row = row - row[pivot] * previous
    positions = row.nonzero_positions()
    if positions:
        pivot = positions[0]
        basis.append((pivot, row / row[pivot]))


# Stage 2
def compute_reductions_dependency(rho0,B,ensure_termination):
    r"""
    Generates and computes a linear dependency among the $\rho_i$ as described in Section 7 of Lairez 2016.

    The matrix $\rho_0$ represents the reconstructed $\rho_0(t)$ with respect to a list of monomials $M$
    in $GF(p)(t)[x_0...x_n,u,v]$, all of which have the form $p(x)u^{n+1}$ (i.e., represent top level forms). 
    The matrix $B$ defines the map $m$ with respect to the basis ``M``.

    This is an internal function for sage_periods, and is not meant to be called by the user.

    INPUT:

    * ``rho0`` -- A Sage matrix over $GF(p)(t)$ representing the reconstructed class $\rho_0(t)$ in the basis $M$.
    * ``B`` -- A square Sage matrix over $GF(p)(t)$ representing the connection action on that basis.
    * ``ensure_termination`` -- Boolean flag to decide whether to use rank specialization filter (no guarantee!)
      or exact rank calculation over $\mathbb F_p(t)$ (guarantee!)

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
    t = K.gen() # Generator for GF(p)(t)
    F = K.base_ring() # GF(p)

    def _draw_point_and_evaluate():
        while True:
            u0 = F.random_element()
            basis = []
            try:
                for v in rho:
                    _extend_echelon(basis, vector(F, [c(u0) for c in v]))
                return u0, basis
            except (ZeroDivisionError, ArithmeticError):
                # Rebuild from scratch if any previous iterate has a pole.
                continue

    m = 0
    rho = [vector(rho0.list())]
    u0, echelon = None, []

    while True:
        verbose("        m: "+str(m),level=1)
        # Note that our degree is guaranteed to stay bounded (i.e. r need not increase)
        # because the M we found in Stage 1 is canonical, as are rho0 and B.
        y = None
        independent = False

        if ensure_termination:
            M = Matrix(K, rho).transpose()
            independent = M.rank() == m + 1
        else:
            try:
                if u0 is None:
                    u0, echelon = _draw_point_and_evaluate()
                else:
                    _extend_echelon(echelon, vector(F, [c(u0) for c in rho[m]]))
            except (ZeroDivisionError, ArithmeticError):
                u0, echelon = _draw_point_and_evaluate()
            independent = len(echelon) == m + 1

        if not independent:
            # Possible dependency; decide exactly over GF(p)(t) by solving for
            # coefficients expressing rho_m in terms of the lower rho_k.
            try:
                rhomat = M if ensure_termination else Matrix(K,rho).transpose()
                y = rhomat[:,:m].solve_right(rho[m])
            except ValueError:
                if not ensure_termination:
                    # Spurious rank drop at u0: rho[m] is independent after all.
                    independent = True
                    u0 = None  # Rebuild all specialized rows on the next iteration.
                else:
                    raise ValueError("Oh no! Couldn't find dependency despite rank drop in exact case.")

        if independent:
            rho.append(vector([ fn.derivative(t)for fn in rho[m]]) + B*rho[m])
            m += 1
            continue

        # At this point, we have our coefficient vector. We could return, but we should clear denominators first to make CRT possible.

        # Clear denominators; return dependency relation (no need to discard primes now)
        Rp = K.ring()
        denoms = [ai.denominator() for ai in y]
        LCM = lcm(denoms) if denoms else Rp.one()
        if not LCM:
            raise ArithmeticError("Dependency denominator LCM unexpectedly vanished.")
        out = []
        for i, (ai, den) in enumerate(zip(y, denoms)):
            quotient, remainder = LCM.quo_rem(den)
            if remainder:
                raise ArithmeticError(
                    f"Dependency denominator {i} does not divide its LCM "
                    f"in characteristic {F.characteristic()}.")
            out.append(ai.numerator() * quotient)
        return out + [LCM]

# Stage 3

# Is now handled by reconstruction.CRTLifter class, directly in compute_period_annihilator.
