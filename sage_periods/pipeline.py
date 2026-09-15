r"""Main stages in our modular pipeline for the reduction algorithm"""

from sage.all import *
from .reconstruction import ReconstructionData, recon_add_rat
from .rham_koszul import RhamKoszulData, gauss_manin_helper

from collections import Counter # For Stage 3

from sage.misc.verbose import verbose, set_verbose



# Core modular loop / state handler
def _modular_reduction_pipeline():

    pass

# Stage 1
def compute_gauss_manin_connection(a,f,r,p):
    r"""
    Compute the Gauss-Manin matrix $B$ for the map defined by closing $\rho_0$ under the map $\rho \mapsto [f^\delta \rho]_r$,
    as described in Section 7.2 of Lairez 2016. Also returns the basis for the space and the projection of $a$ into this space.

    This function is typically called as a subroutine of ``compute_period_annihilator``.

    INPUT:

    * ``a`` -- A polynomial representing the numerator of the function under consideration.
    * ``f`` -- A square-free polynomial representing the square-free denominator of the function under consideration.
    * ``r`` -- Order of the reduction $[\cdot]_r$ which we compute.
    * ``p`` -- The characteristic of the field which we evaluate ``a`` and ``f`` into to do our reductions.
    
    OUTPUT: 
    
    * $M$, $\rho_0(t)$ and $B(t)$ as described in Section 7.2 of Lairez 2016.

    """
    A = a.parent()
    t = A.base_ring().gen()
    bad_points = []

    # Build the K = GF(p)(t) we'll need for reconstruction
    K = PolynomialRing(GF(p),t).fraction_field()
    F = K.base_ring() # GF(p)
    AF = A.change_ring(F) # This should equal RK.A across all RK

    # Reconstruction object for B and rho0
    H = ReconstructionData()

    # Counter for seeing if we must increase r, or if we just had bad evaluation point
    rtoosmall = 0
    # Counter for seeing how many points we've evaluated; useful for breaking out if we suspect a bad prime
    uctr = 0
    gauss_manin_helper_times = []

    # Main loop
    while True:
        u = ZZ.random_element(1, p)   # random evaluation point
        u_QQ = QQ(u)                  # for substitution into coefficients over QQ(t)
        u_K  = K(u)                   # for interpolation / storage mod p
        if u in H.points or u in bad_points:
            continue
        verbose("        u: "+str(u),level=1)

        # Heuristically decide if we've encountered a prime that causes a
        # degenerate specialization. If we end up running a TON of points
        # for an example, then we label this prime as bad and break out.
        uctr += 1

        # Try evaluating f, a, f^delta into our ring. If this fails, pick a new point.
        fdict = f.dict()
        fdelta = A({ key: fdict[key].derivative() for key in fdict.keys()}) #f^delta
        
        try:
            feval = AF(SR(f).subs({t:u_QQ}))
            fdeltaeval = AF(SR(fdelta).subs({t:u_QQ}))
            aeval = AF(SR(a).subs({t:u_QQ}))
        except:
            verbose("Evaluations of f, a, or f^delta failed. Pick a different point.",level=1)
            bad_points.append(u)
            continue
            
        # Build our RhamKoszulData object with evaluated f
        U = RhamKoszulData(feval,r=r)
        
        # Run gauss_manin_helper on <fteval,[aeval]>, with our prime and r
        try:
            ret = gauss_manin_helper(U, fdeltaeval, [aeval], None)
        except RuntimeError as e:
            if str(e) == "INCREASE_R":
                rtoosmall += 1
                bad_points.append(u)
                if rtoosmall >= 3:
                    raise RuntimeError("INCREASE_R")
                continue
            else:
                print("Error in compute_gauss_manin_connection: ")
                print(str(e))
                raise

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

# Stage 2
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


# Stage 3

# TODO: Much of this function body could be delegated to general CRT helper methods
# in reconstruction.py. AI version develops incremental CRT engine, for use with
# polynomial coefficient vectors. Copy --> statically type, and put this portion
# in reconstruction, but keep the "operator" portion here?

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