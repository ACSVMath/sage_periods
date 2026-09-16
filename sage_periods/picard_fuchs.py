r"""
Main functions to compute annihilating D-finite operators for diagonals and periods.
"""

# Sage package imports
from sage.all import *

# Imports from our modules
# from .reconstruction import ReconstructionData, compute_reductions_dependency, lift_operator_across_primes, recon_add_rat
# from .rham_koszul import RhamKoszulData, gauss_manin_helper
from .diagonal import minimize_diagonal_annihilator, diagonal_to_period, normalize_diagonal_arguments
from .preparation import compute_homogenization, compute_prepared_fraction, minimize_denominator_degree, probe_reduction_order
from .pipeline import compute_gauss_manin_connection, compute_reductions_dependency, lift_operator_across_primes
from .errors import SagePeriodsError, ReductionOrderTooSmallError, BadPrimeError, NeedMorePrimesError, ProbeBasisCapExceededError

# Check if ore_algebra is available, and import it if so
from . import _is_ore_algebra_installed
if _is_ore_algebra_installed:
    from ore_algebra import OreAlgebra
else:
    print("Warning: 'ore_algebra' is not installed. Operators will be output as dense univariate OrePolynomials in Q[t][Dt]. Some features may be disabled.")

from sage.misc.verbose import verbose, set_verbose


def compute_diagonal_annihilator(R, r = None, vari = None, Dt = None, t = None, reduction_order = None, minimize = False, certify = False, minimize_denom_deg = None, ncpus = None, ensure_termination = False):
    r"""
    Given a symbolic rational function $R(x_1,...,x_d)$, compute a D-finite equation annihilating the $r$-diagonal of $R$.

    INPUT:

    * ``R`` -- A rational function, either an element of ``SymbolicRing`` or the Fraction Field of a Multivariate Polynomial Ring.
    * ``r`` -- (Optional) A vector of integers specifying the direction of the diagonal, taken as the all ones vector if not specified.
    * ``vari`` -- (Optional) A vector of all variables appearing in ``R``, used to fix the order of the variables for computation. Taken as ``R.variables()`` if not specified.
    * ``Dt`` -- (Optional) The name for the differential symbol in the output, taken as ``Dt`` by default.
    * ``t`` -- (Optional) The name for the variable in the output, taken as ``t`` by default.
    * ``reduction_order`` -- (Optional) The smallest order to start computing the modular Lairez reduction at. When ``None`` is passed in, a small order is selected
     which is likely to give the best combination of small operator size and small runtime.
    * ``minimize`` -- (Optional) When ``True`` (and``ore_algebra`` is installed), 
      recurrence guessing is used to attempt to minimize the operator produced 
      by vanilla ``compute_diagonal_annihilator``. If this function returns with
      ``minimize=True``, this certifies that the output annihilates
      the diagonal and is a right divisor of the operator computed when 
      ``minimize`` is set to ``False``.
    * ``certify`` -- (Optional) When ``True``,
      the partial certificates from Section 7.3 of Lairez are computed and verified
      for operator annihilating the diagonal, and the pair ``(operator, certificate)`` 
      is returned. When combined with ``minimize=True``, the returned operator is 
      certified relative to the pre-minimization operator.
    * ``minimize_denom_deg`` -- (Optional) Boolean flag enabling torus
      change-of-variables preprocessing; see ``compute_period_annihilator``. Can
      greatly speed up computation.
    * ``ncpus`` -- (Optional) Parameter specifying number of CPU cores to use for
      multiprocessing.
    * ``ensure_termination`` -- (Optional) Some speedups (e.g. rank specialization in Stage 2) lack termination guarantees,
      but failure is extremely unlikely. Setting to ``True`` forces execution to fall back to slower methods with guaranteed termination.

    OUTPUT:

    * ``L`` -- An element of the differential ``OreAlgebra`` in ``t`` and ``Dt`` that annihilates the $r$-diagonal of ``R`` or, if ``ore_algebra`` is not available, an element of the ``OrePolynomialRing`` over ``QQ[t]`` with derivation symbol ``Dt``
    * TODO: Add output and EXAMPLES for certificates and other parameters!

    EXAMPLES:

    Computing the main diagonal annihilator on a typical example.

        sage: var('t x y')
        sage: F = 1/(1-x-y-x*y**3)
        sage: compute_diagonal_annihilator(F)
        (t^6 + 2/3*t^5 + 5/9*t^4 + 4/9*t^3 - 1/9*t^2 - 2/81*t + 1/243)*Dt^2 + (4*t^5 + 10/3*t^4 + 8/9*t^2 - 4/27*t + 2/81)*Dt + 2*t^4 + 2*t^3 - 2/3*t^2 - 2/27*t - 8/81

    You can pass nonstandard directions with positive integer entries,

        sage: r = [2,3]
        sage: compute_diagonal_annihilator(F,r=r,vari=[x,y])
        (243*t^8 - 810*t^7 + 837*t^6 - 2412*t^5 - 1779*t^4 - 202*t^3 + 27*t^2)*Dt^3 + (2187*t^7 - 4617*t^6 + 972*t^5 - 594*t^4 - 5661*t^3 - 1557*t^2 + 54*t)*Dt^2 + (4374*t^6 - 3888*t^5 - 4374*t^4 + 3168*t^3 + 2010*t^2 - 1296*t + 6)*Dt + 1458*t^5 + 486*t^4 - 1296*t^3 - 1368*t^2 + 750*t - 30

    and also different names for ``t`` and ``Dt``.

        sage: var('z')
        sage: compute_diagonal_annihilator(F,r=r,vari=[x,y],t=z, Dt='Dz')
        (243*z^8 - 810*z^7 + 837*z^6 - 2412*z^5 - 1779*z^4 - 202*z^3 + 27*z^2)*Dz^3 + (2187*z^7 - 4617*z^6 + 972*z^5 - 594*z^4 - 5661*z^3 - 1557*z^2 + 54*z)*Dz^2 + (4374*z^6 - 3888*z^5 - 4374*z^4 + 3168*z^3 + 2010*z^2 - 1296*z + 6)*Dz + 1458*z^5 + 486*z^4 - 1296*z^3 - 1368*z^2 + 750*z - 30

    If the ``ore_algebra`` package is not installed, a warning will be displayed
    and results will be returned as an ``OrePolynomial`` in ``t`` and ``Dt``

            sage: from sage_periods import compute_diagonal_annihilator
            sage: var('t x y')
            sage: F = 1/(1-x-y-x*y**3)
            sage: L = compute_diagonal_annihilator(F)
            sage: L
            Warning: 'ore_algebra' is not installed. Operators will be output as dense univariate OrePolynomials in Q[t][Dt]. Some features may be disabled.
        (t^6 + 2/3*t^5 + 5/9*t^4 + 4/9*t^3 - 1/9*t^2 - 2/81*t + 1/243)*Dt^2 + (4*t^5 + 10/3*t^4 + 8/9*t^2 - 4/27*t + 2/81)*Dt + 2*t^4 + 2*t^3 - 2/3*t^2 - 2/27*t - 8/81
            sage: L.parent()
            Ore Polynomial Ring in Dt over Univariate Polynomial Ring in t over Rational Field twisted by d/dt

    !!! warning
    
        In particular, ``OrePolynomial`` has no reliable built-in for "evaluating" on elements of the base ring. For instance, calling ``L(t^2 + 2*t - 3)`` will not apply ``L`` to the polynomial ``t^2 + 2*t - 3``, even if both objects are in the right rings. It will instead return an error.
        This is one of several limitations of the ``OrePolynomialRing`` class, hence why we recommend ``OreAlgebra``.

    """

    # Validate inputs and "normalize" arguments so R is symbolic, r and vari are defined
    # (with vari consisting of "standardized" symbolic variables that can be cleared at the end of compute_period_annihilator),
    # vari / r - entry order is optimized for minimal complexity in the resulting diagonal transformation,
    # and t and Dt are defined.
    R_normalized, r, vari_normalized, Dt, t = normalize_diagonal_arguments(R,r,vari,Dt,t)
    
    # !! vari and r may be reordered in order to make our diagonal transformation "nicer."
    # Also, if R was an element in a fraction field, we will have created symbolic variables
    

    # Corner case: If R is constant then either:
    # - We didn't pass in r, in which case we want the main diagonal.
    # - R is its own diagonal, if we passed in r = (1...1)
    # - R's diagonal is zero, if we passed in different r.
    if R_normalized.numerator().is_constant() and R_normalized.denominator().is_constant():
        if r == None or all(r[i] == 1 or r[i] == Integer(1) for i in range(d)):
            return compute_period_annihilator(R_normalized, t, Dt)
        else:
            return compute_period_annihilator(SR(0), t, Dt)
        
    # Compute diagonal transformation.
    G = diagonal_to_period(R_normalized,r,vari_normalized,t)

    # Compute period operator and/or certificates
    L = compute_period_annihilator(G, t, Dt,ensure_termination=ensure_termination)
    # TODO: Add capability for returning certificate.

    # TODO: Add routine for verifying computed certificate.

    # Clear symbolic variables that we created in normalize_diagonal_arguments.
    for v in vari_normalized:
        reset(str(v))

    # TODO: Add possiblity of minimizing the diagonal operator using "guess and verify,"
    # plus rigourous verification of the minimization.
    L_min = minimize_diagonal_annihilator

    return L

    
def compute_period_annihilator(R, t, Dt,reduction_order = None, certify = False, minimize_denom_deg = None, ncpus = None,ensure_termination=False):
    r"""
    Given a symbolic rational function $R(t,x_1,...,x_n)$, compute a D-finite equation 
    annihilating the period integrals (i.e., residues) of $R$ with respect to $x_1,...,x_n$.

    INPUT:

    * ``R`` -- A rational function in $\mathbb{Q}(t)(x_1,...,x_n)$, either an element of ``SymbolicRing`` or the Fraction Field over a Multivariate Polynomial Ring in variables including $t$.
    * ``t`` -- A variable (either symbolic or a generator) appearing in `R` which names the output.
    * ``Dt`` -- A string used to name the operator for differentiation with respect to ``t``.
    * ``reduction_order`` -- (Optional) The smallest order to start computing the modular Lairez reduction at. When ``None`` is passed in, a small order is selected
     which is likely to give the best combination of small operator size and small runtime.
    * ``certify`` -- (Optional) When ``True``,
      the partial certificates from Section 7.3 of Lairez are computed and verified
      for operator annihilating the diagonal, and the pair ``(operator, certificate)`` 
      is returned. When combined with ``minimize=True``, the returned operator is 
      certified relative to the pre-minimization operator.
    * ``minimize_denom_deg`` -- (Optional) Boolean flag enabling torus
      change-of-variables preprocessing; see ``compute_period_annihilator``. Can
      greatly speed up computation.
    * ``ncpus`` -- (Optional) Parameter specifying number of CPU cores to use for
      multiprocessing.
    * ``ensure_termination`` -- (Optional) Some speedups (e.g. rank specialization in Stage 2) lack termination guarantees,
      but failure is extremely unlikely. Setting to ``True`` forces execution to fall back to slower methods with guaranteed termination.

    OUTPUT:

    * ``L`` -- An element of the differential ``OreAlgebra`` in ``t`` and ``Dt`` that 
    annihilates the residue of ``R`` with respect to $x_1,...,x_n$ or, if ``ore_algebra`` 
    is not available, an element of the ``OrePolynomialRing`` over ``QQ[t]`` with derivation symbol ``Dt``
    * TODO: Add output and EXAMPLES for certificates and other parameters!


    EXAMPLES:

    Computing a period annihilator for a typical example.

            sage: var('t x y')
            sage: F = (1+y*x+t)/(1-x*t+x/y)
            sage: compute_period_annihilator(F,t,'Dt')
            (t^4 + t^3 + 3*t)*Dt + t^3 + 2*t^2 + 12

    The result may change significantly depending on which variable is treated as the paramater.

            sage: compute_period_annihilator(F,x,'Dx')
            x*Dx + 1

    Another example, taken from Lairez 2016.

            sage: # Apery example 
            sage: var('w z t')
            sage: F = 1/(1-(1-x*y)*z-t*x*y*z*(1-x)*(1-y)*(1-z))
            sage: compute_period_annihilator(F,t,'Dt')
            (t^8 - 257551/5980*t^7 + 434946/1495*t^6 + 374641/598*t^5 - 528511/1495*t^4 + 58941/5980*t^3)*Dt^4 + (170/13*t^7 - 2993981/5980*t^6 + 5837574/1495*t^5 + 8354793/2990*t^4 - 3626052/1495*t^3 + 54575/1196*t^2)*Dt^3 + (565/13*t^6 - 8477641/5980*t^5 + 68078253/5980*t^4 - 4844873/5980*t^3 - 20736217/5980*t^2 + 6549/230*t)*Dt^2 + (475/13*t^5 - 2841061/2990*t^4 + 10195067/1495*t^3 - 92604/23*t^2 - 15836/23*t - 2183/598)*Dt + 53/13*t^4 - 223701/2990*t^3 + 864791/2990*t^2 - 815887/2990*t + 10915/598

    If the ``ore_algebra`` package is not installed, a warning will be displayed
    and results will be returned as an ``OrePolynomial`` in ``t`` and ``Dt``:

            sage: from sage_periods import compute_period_annihilator
            sage: var('t x y')
            sage: F = 1/(1-x-y-x*y**3)
            sage: L = compute_period_annihilator(F,x,'Dx')
            sage: L
            Warning: 'ore_algebra' is not installed. Operators will be output as dense univariate OrePolynomials in Q[t][Dt]. Some features may be disabled.
            (x^5 - 7/3*x^4 + 5/3*x^3 - 5/27*x^2 - 4/81*x)*Dx^2 + (4*x^4 - 20/3*x^3 + 10/3*x^2 - 20/27*x - 2/81)*Dx + 2*x^3 - 2*x^2 + 2/3*x - 2/27
            sage: L.parent()
            Ore Polynomial Ring in Dx over Univariate Polynomial Ring in x over Rational Field twisted by d/dx

    !!! warning
    
        In particular, ``OrePolynomial`` has no reliable built-in for "evaluating" on elements of the base ring. For instance, calling ``L(t^2 + 2*t - 3)`` will not apply ``L`` to the polynomial ``t^2 + 2*t - 3``, even if both objects are in the right rings. It will instead return an error.
        This is one of several limitations of the ``OrePolynomialRing`` class, hence why we recommend ``OreAlgebra``.

    !!! note

        In the case where $R$ and $t$ are elements of a ``FractionField``, this function will create ``SymbolicRing`` variables corresponding to
        the variables of $R$ and $t$. But these symbolic variables may be shadowed by their fraction field counterparts.

    """
    verbose(f"Computing operator annihilating residue of {R}", level=1)
    
    '''
        TODO:
        * Add criteria for automatically enabling heuristics for variable order etc., based on the complexity of R.
        * Build fraction field with t in coefficient field *early.* That way we can use it later on and reap the speedups.
        * Collect all algebraic irrational coefficients found in R and build a finite algebraic extension over QQ, 
        then use this to make K.
    '''

    # If we've passed in a non-SR R, verify assumptions on it, then put into SR.
    if R.parent() != SR:
        assert t in R.parent().gens(), "t must be a generator of the parent ring of R."
        t = SR.var(t)
        R = SR(R)
    else:
        if R.parent() == SR:
            assert t.parent() == SR, "If R is symbolic, then t must be a symbolic ring element as well."
    vari = [v for v in R.variables() if v != t]
    # # Normalize space variables -- just in case they're not from a diagonal.
    # vari_normalized = list(SR.var("sage_periods_x", len(vari)))
    # R_normalized = R.subs({old:new for (old,new) in zip(vari,vari_normalized)}) # This might not be necessary, here.
    
    # Define base field and set up variables
    F = QQ  # TODO: Add Possiblity of working over algebraically defined coefficients.
    K = PolynomialRing(F,t).fraction_field()
    t_symbolic = t
    t = K.gen()

    # Introduce homogenizing variable
    extra_var = SR.var('sage_periods_h')
    A = PolynomialRing(K,vari+[extra_var],len(vari)+1, order="degrevlex")
    B = A.fraction_field()

    # Cast R into B; free symbolic variables.
    R = B(R)
    # reset(str(t_symbolic))
    # for x in vari_normalized:
    #     reset(str(x))
    
    # Prepare our rational function
    Fhom = compute_homogenization(R)
    a,f,q = compute_prepared_fraction(Fhom)
    verbose(f"The prepared fraction has the form (a,f,q) = {(a,f,q)}",level=1)


    # TODO: Run degree minimization on the denominator, if it is called


    # Build OreAlgebra object associated with this ring (if possible)
    # First, recast A as the poly ring in t over F[x_0,...,x_n]
    _A_iso = PolynomialRing(F,t)
    _t_iso = _A_iso.gen()
    
    # Now build the correct structure with t as the variable and Dt acting on it
    if _is_ore_algebra_installed:
        Alg = OreAlgebra(_A_iso, (Dt, lambda p: p,
                                    lambda p: p.derivative(_t_iso))) # This allows us more customization in the naming of our operator variable.
        Dt = Alg.gen()
    else:
        Alg = PolynomialRing(_A_iso, Dt,order='degrevlex')
        Alg = OrePolynomialRing(_A_iso, _A_iso.derivation(), Dt)
        Dt = Alg.gen()

    verbose("Starting to compute Picard-Fuchs operator using evaluation-interpolation.",level=1)

    # Choose the reduction order.  If the user passed one, start there.
    # If none was passed, probe r = 1 and r = 2 modulo one prime
    # and start at whichever produces the operator of smaller order, 
    # if both execution times are under 30 seconds.
    # Else, pick whichever one is faster.
    seed = {}
    profile = {}
    if reduction_order is None:
        verbose("Choosing the reduction order by probing modulo one prime.",level=1)
        # print(probe_reduction_order(a,f))
        r, seed = probe_reduction_order(a, f,ensure_termination) # Future: add , "ncpus=ncpus", passed in from caller. Also, probe_reduction_order should return "profile" as well.
        verbose(f"Auto-probe selected reduction order r = {r}.",level=1)
    else:
        r = Integer(reduction_order)
        assert r >= 1, "reduction_order must be a positive integer."

    # Run a modular algorithm to compute L mod p. If no relations are found, r will be increased.
    # If the denominator of R defines a smooth variety, the algorithm is guaranteed to terminate with r=1. In general the algorithm
    # will terminate, but the smallest value of r guaranteed to make it stop is still an open problem.
    deq = None
    while not deq:
        verbose("    r: " + str(r),level=1)
        try:
            ### Run the modular algorithm, with fixed r ###
            # Hold coefficients for our Picard-Fuchs operator over F_p(t) in a dictionary keyed by p
            Lp_coeffs_dict = {}
            # Store primes p which can't be used in CRT lifting
            bad_primes = [] 
            primectr = 0

            # Buckets for reconstructed operators (modulo many primes)
            L_recon_candidate = None
            '''
                Notes on picking primes:
                - Should be less than 2^31, so that we can store field elements as single int64 objects
                - Should be less than 2^29, because of this error:
                ``NotImplementedError: Division of multivariate polynomials over prime fields with characteristic > 2^29 is not implemented.``
                (since we divide by RK.w_prime)

                - According to Lairez 2016, we should pick a prime much larger than e*n*N to give a high probability of not yielding a degenerate specialization.
                - Here we take a lower bound that's a couple of orders of magnitude lower than our upper bound.
            '''
            p = random_prime(2**26 - 1, proof=True, lbound=2**23) # There are over 3 million primes in this interval
            while True:
                # Find a prime that's not in our list yet
                while p in Lp_coeffs_dict or p in bad_primes:
                    p = random_prime(2**26-1, proof=True, lbound=2**23)

                verbose("p: "+str(p),level=1)
                primectr +=1
        
                verbose(f"Computing differential operator modulo prime number {primectr} ({p}).",level=1)
        
                #######################################################################################
                # STAGE 1: Compute M, rho_0(t) and B(t) through evaluation and rational interpolation. 
                # This is the only stage where it might be necessary to increase r.
                #######################################################################################
                try:
                    rho0, M, B = compute_gauss_manin_connection(a,f,r,p)
                except ReductionOrderTooSmallError as e:
                    # Increase r (we're inside a try already)
                    raise e
                
                # If we get a trivial basis, then we're done.
                if len(M) == 0 or rho0.nrows() == 0:
                    return Alg(1)
                    
                ###############################################################################################
                # STAGE 2: Use matrix formula for rho[i] to compute a dependency among the rho[i].
                # Note that compute_reductions_dependency returns ALL coeffs up to m, with cleared denominators.
                # (in particular this means rho[m] won't usually have 1 as its leading coeff, so our operator isn't monic.)
                ###############################################################################################
                verbose("Computing linear relation... ",level=1)
                try:
                    Lp_coeffs_denoms_cleared = compute_reductions_dependency(rho0,B,ensure_termination=ensure_termination)
                except BadPrimeError:
                    bad_primes.append(p)
                    continue
                
                Lp_coeffs_dict[p] = Lp_coeffs_denoms_cleared
                
                # At this point, we have a list of coeffs for our forms rho[i], evaluated in F_p(t), constituting an operator L_p over F_p(t)
                ###############################################################################################################
                # STAGE 3:  Use Chinese remainder theorem to lift these coeffs a_i to QQ(t).
                # The list of operators L_p is encoded by Lp_coeffs_dict.
                # Assumes we cleared denominators across lists at end of Stage 2 so that Lp_coeffs_dict contains polynomials.
                ###############################################################################################################
                verbose(f"Found an equation of order {len(Lp_coeffs_denoms_cleared)-1} and degree {max([p.numerator().degree() for p in Lp_coeffs_denoms_cleared])}.",level=1)
                try:
                    L_coeffs = lift_operator_across_primes(t,Lp_coeffs_dict,bad_primes)
                except NeedMorePrimesError:
                    # Either we pruned out too many primes in CRT, or we don't have enough to rationally reconstruct. Add more.
                    verbose("Rational reconstruction failed to lift one or more of our coefficients. Need more primes.",level=1)
                    continue

                # Check to see if L_coeffs has been assigned in L_recon_buckets and if it equals this value (stability)
                if L_recon_candidate is not None and L_recon_candidate == L_coeffs:
                    # We take this as the reconstructed operator having stabilized
                    # Break and return the operator
                    break
                else:
                    L_recon_candidate = L_coeffs
                    continue
            
        
            # Build operator then return it.
            m = len(L_coeffs)
            deq = Alg(L_coeffs[m-1]*Dt**(m-1) - sum([L_coeffs[k]*(Dt**k) for k in range(m-1)]))
        except ReductionOrderTooSmallError as e:
            r += 1
            continue

    return deq
