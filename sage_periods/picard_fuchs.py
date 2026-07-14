r"""
Main functions to compute annihilating D-finite operators for diagonals and periods.
"""

# Sage package imports
from sage.all import *
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.ideal import Ideal
from sage.rings.integer import Integer
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.symbolic.ring import SR
from sage.rings.polynomial.ore_polynomial_ring import OrePolynomialRing
from sage.misc.reset import reset
from sage.rings.integer_ring import ZZ
from sage.misc.misc_c import prod
from sage.arith.misc import random_prime
from sage.calculus.var import var
from sage.sets.set import Set

# Imports from our modules
from .helpers import first_coordinate_section
from .reconstruction import ReconstructionData, compute_reductions_dependency, lift_operator_across_primes, recon_add_rat
from .rham_koszul import RhamKoszulData, gauss_manin_helper

# Check if ore_algebra is available, and import it if so
from . import _is_ore_algebra_installed
if _is_ore_algebra_installed:
    from ore_algebra import OreAlgebra
else:
    print("Warning: 'ore_algebra' is not installed. Operators will be output as dense univariate OrePolynomials in Q[t][Dt]. Some features may be disabled.")

from sage.misc.verbose import verbose, set_verbose


def _prepared_fraction_data(R):
    """Return ``(a, f, q, factors)`` for a fraction-field element ``R``."""
    A = R.parent().ring()
    if R.numerator().degree() <= 0 and R.denominator().degree() == 0:
        return A(R), A.one(), Integer(1), []

    factorization = list(R.denominator().factor())
    q = max(exponent for _, exponent in factorization)
    factors = [factor for factor, _ in factorization]
    f = A(prod(factors))
    return A(f**q * R), f, q, factors


def _diagonal_integrand(R, r, vari, pivot_index, t):
    """Construct the diagonal period integrand for one pivot."""
    order = [pivot_index] + [i for i in range(len(r)) if i != pivot_index]
    ordered_vari = [vari[i] for i in order]
    ordered_r = [r[i] for i in order]
    pivot, *nonpivot_vari = ordered_vari
    pivot_weight = ordered_r[0]

    if pivot_weight == 1:
        weighted_product = prod(
            x**ri for x, ri in zip(nonpivot_vari, ordered_r[1:])
        )
        G = R.subs({pivot: pivot / weighted_product})
    else:
        G = first_coordinate_section(
            R, pivot, pivot_weight, u=pivot, vari=ordered_vari
        )
        weighted_product = prod(
            x**(pivot_weight * ri)
            for x, ri in zip(nonpivot_vari, ordered_r[1:])
        )
        G = G.subs({
            pivot: pivot / weighted_product,
            **{x: x**pivot_weight for x in nonpivot_vari},
        })

    G /= prod(nonpivot_vari)
    return (
        G.subs({pivot: t, t: pivot})
        if t in ordered_vari and t != pivot
        else G.subs({pivot: t})
    )


def _pivot_candidate(G, pivot, t, K, homogenizing_variable):
    """Prepare one pivot and compute its static complexity metrics."""
    integration_variables = sorted(
        (set(G.numerator().variables()) | set(G.denominator().variables())) - {t},
        key=str,
    )
    if homogenizing_variable in integration_variables:
        raise ValueError(
            "Please pick a different variable name. Can't use extra_var in "
            "input R, since it's reserved for homogenization."
        )

    names = integration_variables + [homogenizing_variable]
    A = PolynomialRing(K, names, len(names), order="degrevlex")
    Fhom = compute_homogenization(A.fraction_field()(G))
    a, f, q, factors = _prepared_fraction_data(Fhom)
    N = f.degree()
    n = A.ngens() - 1
    fdict = f.dict()
    fdelta = A({
        monomial: derivative
        for monomial, coefficient in fdict.items()
        if (derivative := coefficient.derivative()) != 0
    })

    leading_coefficient = f.lc()
    normalized_coefficients = [
        coefficient / leading_coefficient for coefficient in fdict.values()
    ]
    t_degree = max(
        max(
            coefficient.numerator().degree(),
            coefficient.denominator().degree(),
        )
        for coefficient in normalized_coefficients
    )
    moving_components = 0
    for factor in factors:
        factor_lc = factor.lc()
        moving_components += any(
            (coefficient / factor_lc).derivative() != 0
            for coefficient in factor.dict().values()
        )

    return {
        "pivot": pivot,
        "integrand": G,
        "a": a,
        "f": f,
        "fdelta": fdelta,
        "static": (
            binomial(q * N - 1, n),
            q,
            N,
            len(fdelta.dict()),
            t_degree,
            moving_components,
        ),
    }


def _pareto_front(candidates):
    """Discard candidates dominated componentwise in their static metrics."""
    def dominates(left, right):
        comparisons = [
            x <= y for x, y in zip(left["static"], right["static"])
        ]
        return all(comparisons) and left["static"] != right["static"]

    return [
        candidate
        for candidate in candidates
        if not any(
            other is not candidate and dominates(other, candidate)
            for other in candidates
        )
    ]


def _evaluate_polynomial_mod_p(poly, target_ring, point, field):
    """Evaluate the coefficient field parameter of ``poly`` modulo ``p``."""
    return target_ring({
        monomial: (
            field(coefficient.numerator()(point))
            / field(coefficient.denominator()(point))
        )
        for monomial, coefficient in poly.dict().items()
    })


def _gauss_manin_pilot(candidate, prime, attempts=12):
    """Return a two-specialization score ``(|M|, nnz(B), ...)`` or ``None``."""
    A = candidate["a"].parent()
    field = GF(prime)
    target_ring = A.change_ring(field)
    results = []
    used_points = set()

    for _ in range(attempts):
        if len(results) == 2:
            break
        point = ZZ.random_element(1, prime)
        if point in used_points:
            continue
        used_points.add(point)

        try:
            feval = _evaluate_polynomial_mod_p(
                candidate["f"], target_ring, point, field
            )
            fdeltaeval = _evaluate_polynomial_mod_p(
                candidate["fdelta"], target_ring, point, field
            )
            aeval = _evaluate_polynomial_mod_p(
                candidate["a"], target_ring, point, field
            )
            result = gauss_manin_helper(
                RhamKoszulData(feval, r=1), fdeltaeval, [aeval]
            )
        except RuntimeError as error:
            if str(error) == "INCREASE_R":
                return None
            raise
        except Exception as error:
            verbose(
                f"Discarding pivot pilot point {point}: {error}", level=2
            )
            continue

        results.append((len(result.basis), len(result.gm.dict()), result.ebasis))

    if len(results) < 2:
        return None

    basis_sizes = [size for size, _, _ in results]
    nonzero_counts = [nonzero for _, nonzero, _ in results]
    return (
        max(basis_sizes),
        max(nonzero_counts),
        int(results[0][2] != results[1][2]),
        sum(basis_sizes),
        sum(nonzero_counts),
    )


def _select_pivot_candidate(candidates, original_pivot, vari):
    """Apply the static Pareto filter and, if needed, the modular pilot."""
    for candidate in candidates:
        verbose(
            f"Pivot candidate {vari[candidate['pivot']]} "
            f"(index {candidate['pivot']}) has static metrics "
            f"{candidate['static']}.",
            level=1,
        )

    finalists = _pareto_front(candidates)
    if len(finalists) == 1:
        return finalists[0]

    pilot_prime = random_prime(2**20 - 1, proof=True, lbound=2**19)
    successful = []
    for candidate in finalists:
        candidate["pilot"] = _gauss_manin_pilot(candidate, pilot_prime)
        verbose(
            f"Pivot candidate {vari[candidate['pivot']]} "
            f"(index {candidate['pivot']}) has pilot score "
            f"{candidate['pilot']}.",
            level=1,
        )
        if candidate["pilot"] is not None:
            successful.append(candidate)

    if not successful:
        return next(
            candidate for candidate in candidates
            if candidate["pivot"] == original_pivot
        )

    best_score = min(candidate["pilot"] for candidate in successful)
    tied = [candidate for candidate in successful if candidate["pilot"] == best_score]
    return next(
        (
            candidate for candidate in tied
            if candidate["pivot"] == original_pivot
        ),
        tied[0],
    )


def compute_diagonal_annihilator(R, r=None, vari=None, Dt=None, t=None):
    r"""
    Given a symbolic rational function $R(x_1,...,x_d)$, compute a D-finite equation annihilating the $r$-diagonal of $R$.

    INPUT:

    * ``R`` -- A rational function, either an element of ``SymbolicRing`` or the Fraction Field of a Multivariate Polynomial Ring.
    * ``r`` -- (Optional) A vector of integers specifying the direction of the diagonal, taken as the all ones vector if not specified.
    * ``vari`` -- (Optional) An ordered vector of all variables appearing in ``R``. It is required for a non-main direction and otherwise defaults to the variables detected in ``R``.
    * ``Dt`` -- (Optional) The name for the differential symbol in the output, taken as ``Dt`` by default.
    * ``t`` -- (Optional) The name for the variable in the output, taken as ``t`` by default.

    OUTPUT:

    * ``L`` -- An element of the differential ``OreAlgebra`` in ``t`` and ``Dt`` that annihilates the $r$-diagonal of ``R`` or, if ``ore_algebra`` is not available, an element of the ``OrePolynomialRing`` over ``QQ[t]`` with derivation symbol ``Dt``

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
    assert Dt != "D", "Please pick a different operator symbol; D is reserved."
    assert Dt is None or isinstance(Dt, str), (
        "Derivative symbol Dt must be None or a string."
    )
    assert not (Dt is not None and t is None), (
        "If you pass Dt, then you must pass t as well."
    )

    if t is None:
        t, Dt = var("t"), "Dt"
    elif Dt is None:
        Dt = f"D{t}"

    input_is_symbolic = R.parent() == SR
    if input_is_symbolic:
        variables_in_R = list(R.variables())
    else:
        used = set(R.numerator().variables()) | set(R.denominator().variables())
        variables_in_R = [x for x in R.parent().gens() if x in used]

    if r is None:
        vari = variables_in_R if vari is None else list(vari)
        r = [Integer(1)] * len(vari)
        is_main_diagonal = True
    else:
        assert len(r) > 0, "Direction vector r must be nonempty."
        assert all(isinstance(ri, (int, Integer)) for ri in r), (
            "Direction vector r must be a list of integers."
        )
        assert all(ri > 0 for ri in r), (
            "Direction vector r must have strictly positive entries."
        )
        r = [Integer(ri) for ri in r]
        is_main_diagonal = all(ri == 1 for ri in r)
        if not is_main_diagonal and gcd(r) > 1:
            print(
                "Warning: r has nontrivial gcd. Ensure you meant to do this, "
                "as computations will be slower and operators larger."
            )
        if vari is None:
            assert is_main_diagonal, (
                "vari must be provided when r is not the all-ones vector."
            )
            vari = variables_in_R
        else:
            vari = list(vari)

    assert len(vari) == len(r), (
        "Direction vector r must have the same length as vari."
    )
    assert (
        len(vari) == len(variables_in_R)
        and len(set(vari)) == len(vari)
        and set(vari) == set(variables_in_R)
    ), "vari must contain each variable appearing in R exactly once."
    assert (
        all(x.parent() == SR for x in vari)
        if input_is_symbolic
        else all(x in R.parent().gens() for x in vari)
    ), "vari contains variables incompatible with R."

    temporary_symbols = []
    if not input_is_symbolic:
        symbolic_generators = {x: var(str(x)) for x in variables_in_R}
        temporary_symbols = list(symbolic_generators.values())
        R = SR(R)
        vari = [symbolic_generators[x] for x in vari]

    try:
        if R.numerator().is_constant() and R.denominator().is_constant():
            return compute_period_annihilator(
                R if is_main_diagonal else SR(0), t, Dt
            )

        minimum_weight = min(r)
        candidate_pivots = [
            i for i, ri in enumerate(r) if ri == minimum_weight
        ]
        original_pivot = candidate_pivots[0]

        if len(candidate_pivots) == 1:
            best_pivot = original_pivot
            best_integrand = _diagonal_integrand(R, r, vari, best_pivot, t)
            verbose(
                f"Using the only eligible pivot {vari[best_pivot]} "
                f"(index {best_pivot}).",
                level=1,
            )
        else:
            K = PolynomialRing(QQ, t).fraction_field()
            homogenizing_variable = var("extra_var")
            candidates = [
                _pivot_candidate(
                    _diagonal_integrand(R, r, vari, pivot, t),
                    pivot,
                    t,
                    K,
                    homogenizing_variable,
                )
                for pivot in candidate_pivots
            ]
            winner = _select_pivot_candidate(
                candidates, original_pivot, vari
            )
            best_pivot = winner["pivot"]
            best_integrand = winner["integrand"]
            verbose(
                f"Using pivot {vari[best_pivot]} (index {best_pivot}).",
                level=1,
            )

        return compute_period_annihilator(best_integrand, t, Dt)
    finally:
        for symbol in temporary_symbols:
            reset(str(symbol))


def compute_period_annihilator(R, t, Dt):
    r"""
    Given a symbolic rational function $R(t,x_1,...,x_n)$, compute a D-finite equation 
    annihilating the period integrals (i.e., residues) of $R$ with respect to $x_1,...,x_n$.

    INPUT:

    * ``R`` -- A rational function in $\mathbb{Q}(t)(x_1,...,x_n)$, either an element of ``SymbolicRing`` or the Fraction Field over a Multivariate Polynomial Ring in variables including $t$.
    * ``t`` -- A variable (either symbolic or a generator) appearing in `R` which names the output.
    * ``Dt`` -- A string used to name the operator for differentiation with respect to ``t``.

    OUTPUT:

    * ``L`` -- An element of the differential ``OreAlgebra`` in ``t`` and ``Dt`` that 
    annihilates the residue of ``R`` with respect to $x_1,...,x_n$ or, if ``ore_algebra`` 
    is not available, an element of the ``OrePolynomialRing`` over ``QQ[t]`` with derivation symbol ``Dt``

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

    """
    verbose(f"Computing operator annihilating residue of {R}", level=1)
    
    '''
        TODO:
        * Add criteria for automatically enabling heuristics for variable order etc., based on the complexity of R.
        * Build fraction field with t in coefficient field *early.* That way we can use it later on and reap the speedups.
        * Collect all algebraic irrational coefficients found in R and build a finite algebraic extension over QQ, 
        then use this to make K.
    '''
    
    # If we're passed in a non-SR R, verify assumptions on it, then put us into the symbolic ring.
    # Laziness
    if R.parent() != SR:
        assert t in R.parent().gens(), "t must be a generator of the parent of R."
        sym_vari = [var(sym) for sym in [ str(x) for x in R.parent().gens() if x != t]]
        t = var(str(t))
        R = SR(R)
    
    # Define base field and set up variables in Sage
    F = QQ  
    K = PolynomialRing(F,t).fraction_field()
    t_symbolic = t
    t = K.gen()
    vari = sorted((set(R.numerator().variables()) | set(R.denominator().variables())) - {t_symbolic},key=str)

    # Introduce homogenizing variable
    extra_var = var('extra_var')
    if extra_var in vari:
        raise Exception("Please pick a different variable name. Can't use extra_var in input R, since it's reserved for homogenization.")
    A = PolynomialRing(K,vari+[extra_var],len(vari)+1, order="degrevlex")
    B = A.fraction_field()

    # Cast R into B; free symbolic variables, if we created them.
    if R.parent() != SR:
        reset(str(t_symbolic))
        for x in sym_vari:
            reset(str(x))
        R = B(R)
    else:
        R = B(R)
    
    # Prepare our rational function
    Fhom = compute_homogenization(R)
    a,f,q = compute_prepared_fraction(Fhom)
    verbose(f"The prepared fraction has the form (a,f,q) = {(a,f,q)}",level=1)

    # Build OreAlgebra object associated with this ring (if possible)
    # First, "recast" A as the poly ring in t over F[x_0,...,x_n]
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

    # Start with parameter r=1 and run a modular algorithm to compute L mod p. If no relations are found, r will be increased.
    # If the denominator of R defines a smooth variety, the algorithm is guaranteed to terminate with r=1. In general the algorithm
    # will terminate, but the smallest value of r guaranteed to make it stop is still an open problem.
    r = 1
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
                except Exception as e:
                    if str(e) == "INCREASE_R":
                        raise Exception("INCREASE_R")
                    else:
                        print("Error in Stage 1 (Computing M, rho_0 and B): ")
                        print(str(e))
                        raise
                
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
                    Lp_coeffs_denoms_cleared = compute_reductions_dependency(rho0,B)
                except Exception as e:
                    if str(e) == "BAD_PRIME":
                        # This means we couldn't clear denoms successfully in Stage 2. Discard this prime and continue
                        bad_primes.append(p)
                        print("Encountered bad prime in Stage 2. p = "+str(p))
                        continue
                    else:
                        print("Error in Stage 2 (Computing a dependency among the rho[i]): ")
                        print(str(e))
                        raise
                
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
                except Exception:
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
        except Exception as e:
            # If our r is too low
            if str(e) == "INCREASE_R":
                r += 1
                continue
            else:
                print("Error in main r-loop:")
                print(str(e))
                raise
    return deq
    

def compute_homogenization(R, k=None, t = None, homog_var = None):
    r"""
    Return the homogenization of a rational function.

    INPUT:

    * ``R`` -- An element of the SymbolicRing or an element of $F(t)(x_1,...,x_n)$ where
        $F$ is an algebraic extension of $\mathbb{Q}$.
    * ``k`` -- (Optional) An integer specifying the degree in which to homogenize ``R``. If not
        provided, ``k`` is taken to be $-n-1$, where $n$ is the
        number of variables appearing in ``R`` other than the parameter ``t``.
    * ``t`` -- (Optional) An element of the ``SymbolicRing`` defining the parameter when ``R`` is symbolic.
    * ``homog_var`` -- (Optional) If provided, this is the variable used to homogenize.
        If ``R`` is symbolic, this variable should not appear in ``R``. Unused if ``R`` is an element of $F(t)(x_1,...,x_n)$,
        where homogenization occurs with respect to the last generator of the ring.

    ASSUMPTIONS:

    * ``t`` is not ``homog_var`` if ``R`` is symbolic.
    * Homogenizing variable is the *last* variable to appear in ``R.parent().gens()`` if ``R`` is not symbolic. 

    OUTPUT:

    * The homogenization of ``R(t)`` in degree ``k`` with homogenizing variable ``homog_var``.
      Concretely, if $R = R(x_1, ..., x_n, t)$ and ``homog_var = h``, then the result is 
      $h^k R(x_1/h, ..., x_n/h)$. If ``R`` is symbolic, then so is the output.

    EXAMPLES:

    If the input *is* symbolic, you must pass in the parameter and homogenization variable.

            sage: var('t x y')
            sage: F = 1/(1-t*x*y + y**2)
            sage: compute_homogenization(F)
            Traceback (most recent call last):
            ...
            AssertionError: You must provide t and homog_var if R is in the Symbolic Ring.

    By default, the homogenization has degree $-n-1$

            sage: var('h')
            sage: compute_homogenization(F,t=t,homog_var=h)
            -1/(h^3*(t*x*y/h^2 - y^2/h^2 - 1))

    but you can change this by passing a value for ``k``.

            sage: compute_homogenization(F,k=5,t=t,homog_var=h)
            -h^5/(t*x*y/h^2 - y^2/h^2 - 1)
            sage: compute_homogenization(F,k=-5,t=t,homog_var=h)
            -1/(h^5*(t*x*y/h^2 - y^2/h^2 - 1))
            
    You can also pass in an element of $\mathbb{Q}(t)(x_0,...,x_n)$. Here, the parameter is extracted as the 
    generator of the base ring and the homogenization variable is $x_n$.

            sage: K = QQ[t].fraction_field()
            sage: A = K[x,y,h].fraction_field()
            sage: t = K.gen()
            sage: x = A.gens()[0]
            sage: y = A.gens()[1]
            sage: F_frac = 1/(1-t*x*y + y**2)
            sage: compute_homogenization(F_frac)
            1/((-t)*x*y*h + y^2*h + h^3)

    In the non-symbolic case you should not pass in values for ``t`` or ``homog_var``.

            sage: compute_homogenization(F_frac,t=t,homog_var = A.gens()[2])
            Traceback (most recent call last):
            ...
            AssertionError: If R is a member of a fraction field, put the homogenizing variable as last generator of R's parent ring.

    """
    assert R.parent() != SR or (t != None and homog_var != None), "You must provide t and homog_var if R is in the Symbolic Ring."
    assert R.parent() == SR or (t == None and homog_var == None), "Only provide t and homog_var if R is in the Symbolic Ring."
    assert not(R.parent() != SR and homog_var != None), "If R is a member of a fraction field, put the homogenizing variable as last generator of R's parent ring."
    assert not(R.parent() == SR and homog_var in R.variables()), "If R is symbolic, then homogenizing variable can't be apart of R already."
    if R.parent() == SR:
        h = homog_var
        vari = [v for v in R.variables() if v != h if v!= t]

    else:
        B = R.parent()
        h = B.gens()[-1]
        h_poly = B.ring().gens()[-1]
        assert R.numerator().degree(h_poly) <= 0 and R.denominator().degree(h_poly) == 0, "Not allowed to have homogenization variable inside your function already."
        vari = [v for v in B.gens() if v != h]
    if( k == None):
        k = -len(vari)-1
    return h**k*R.subs({ x : x/h for x in vari})


def compute_prepared_fraction(R,t = None):
    r"""
    Given a rational function $R$, either symbolic or in $F(t)(x_0,...,x_n)$ 
    with $F$ a finite algebraic extension of $\mathbb{Q}$, return polynomials $a$ 
    and $f$ and integer $q$ such that $R = a/f^q$. 

    Note that $f$ and $a$ may not be coprime, but $f$ is squarefree.

    INPUT:

    * ``R`` -- An element of the symbolic ring or $F(t)(x_0,...,x_n)$ with $F$ 
        a finite algebraic extension of $\mathbb{Q}$.
    *   ``t`` -- (Optional) A symbolic variable which is the parameter of ``R``, when ``R`` is symbolic.

    OUTPUT:

    * A triple ``(a, f, q)`` such that ``R = a/f^q`` and ``f`` is
      squarefree. If ``R`` is homogeneous then both ``a`` and ``f``
      are homogeneous polynomials. If R is symbolic, then so are ``a`` and ``f``.

    EXAMPLES:

    If the input is symbolic, you must pass in the parameter ``t``.

            sage: var('t x y')
            sage: F = 1/(1-t*x*y + y**2)
            sage: compute_prepared_fraction(F)
            Traceback (most recent call last):
            ...
            AssertionError: If R is symbolic, please provide the parameter t.

            sage: compute_prepared_fraction(F,t)
            (-1, t*x*y - y^2 - 1, 1)

    Like ``compute_homogenization``, the user is allowed to pass in an element
    of a fraction field.

            sage: var('h')
            sage: K = QQ[t].fraction_field()
            sage: A = K[x,y,h].fraction_field()
            sage: t = K.gen()
            sage: x = A.gens()[0]
            sage: y = A.gens()[1]
            sage: F_frac = 1/(1-t*x*y + y**2)
            sage: compute_prepared_fraction(F_frac)
            (-1, t*x*y - y^2 - 1, 1)

    """
    is_symbolic = R.parent() == SR
    assert not (not is_symbolic and t is not None), (
        "Only provide t if R is in the Symbolic Ring."
    )
    assert not (is_symbolic and t is None), (
        "If R is symbolic, please provide the parameter t."
    )

    if is_symbolic:
        t_symbolic = t
        K = PolynomialRing(QQ, t).fraction_field()
        variables = sorted(
            (
                set(R.numerator().variables())
                | set(R.denominator().variables())
            ) - {t_symbolic},
            key=str,
        )
        A = PolynomialRing(K, variables, len(variables), order="degrevlex")
        R = A.fraction_field()(R)

    a, f, q, _ = _prepared_fraction_data(R)
    return (SR(a), SR(f), q) if is_symbolic else (a, f, q)


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

    fdict = f.dict()
    fdelta = A({ key: fdict[key].derivative() for key in fdict.keys()}) #f^delta

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


# __all__ = [
#     "compute_homogenization",
#     "compute_prepared_fraction",
#     "compute_period_annihilator",
#     "compute_diagonal_annihilator",
#     "compute_gauss_manin_connection"
# ]