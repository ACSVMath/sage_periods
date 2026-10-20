r"""Functions used in preparing rational period intetgral for modular pipeline."""

from sage.all import *

from .pipeline import compute_gauss_manin_connection, compute_reductions_dependency, prepare_coefficient_data
from .errors import *
from .profiling import ProfileData
from .reconstruction import reconstruction_degree_hint
from .rham_koszul import _staircase_cache

from sage.misc.verbose import verbose
import multiprocessing
import traceback
import random
import time

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
    if R.parent() == SR:
        is_symbolic = True
    else:
        is_symbolic = False
    
    assert not (not is_symbolic and t != None), "Only provide t if R is in the Symbolic Ring."
    assert (not is_symbolic or t != None), "If R is symbolic, please provide the parameter t."
    if is_symbolic:
        # Need to build the fraction ring and evaluate R
        F = QQ 
        K = PolynomialRing(F,t).fraction_field()
        t_symbolic = t
        t = K.gen()
        vari = sorted((set(R.numerator().variables()) | set(R.denominator().variables())) - {t_symbolic},key=str)
        A = PolynomialRing(K,vari,len(vari), order="degrevlex")
        B = A.fraction_field()
        R = B(R)
    else:
        B = R.parent()
        A = B.ring()

    # Corner case: Constant R. Need this so max() doesn't throw an error.
    if R.numerator().degree() <= 0 and R.denominator().degree() == 0:
        a = A(R)
        f = A(1)
        q = 1
        if is_symbolic:
            return SR(a), SR(f), q
        else:
            return a,f,q

    vari = B.gens()
    g = R.numerator()
    h = R.denominator()

    # Make f squarefree
    h_factors = h.factor()
    pairs = list(h_factors)
    unit = h_factors.unit()
    
    q = max(e for _, e in pairs)
    f = prod(p for p, _ in pairs)
    a = f**q * R 
    if is_symbolic:
        return SR(a), SR(f), q
    else:
        return A(a),A(f),q

def minimize_denominator_degree(R, xgens, rounds=3, concurrent=10, time_cap=30.0):
    r"""
    Torus (monomial) change of variables lowering the denominator degree.

    A random walk over the unimodular substitutions $x_i \mapsto 1/x_i$ and
    $x_i \mapsto x_i x_j^{\pm 1}$ ($x_i, x_j$ in ``xgens``), applied to the
    $dx/x$-normalized integrand $h = R\,\prod_{x_i \in xgens} x_i$, keeping a
    pool of the best candidates and stopping after ``rounds`` rounds without
    improvement (or after ``time_cap`` seconds).  Returns the transformed
    integrand $h'/\prod_{x_i \in xgens} x_i$, or ``R`` itself if no
    improvement was found.

    ``concurrent'' measures how many items we should keep in our pool concurrently.

    This is an internal function for sage_periods, and is not meant to be
    called by the user.
    """
    rng = random.Random(20260803)
    if R.is_zero():
        return R
    B = R.parent()
    A = B.ring()
    h = R * prod(xgens)
    num_coeffs = list(h.numerator().dict().items())
    den_coeffs = list(h.denominator().dict().items())
    num_exps = tuple(tuple(e) for e, _ in num_coeffs)
    den_exps = tuple(tuple(e) for e, _ in den_coeffs)

    # Indicator vector of the torus coordinates divided out by prod(xgens).
    gen_index = {str(g): k for k, g in enumerate(A.gens())}
    xidx = {gen_index[str(v)] for v in xgens}
    torus = tuple(1 if k in xidx else 0 for k in range(A.ngens()))

    # Restrict the walk to the variables that actually occur in h: the
    # integrand may be independent of some torus coordinates.A monomial 
    # substitution in an absent variable fixes h, so we needn't include them.
    walk = [k for k in sorted(xidx)
            if any(e[k] for e in num_exps) or any(e[k] for e in den_exps)]
    if not walk:
        return R

    st = (num_exps, den_exps)
    start = _denominator_degree(st, torus)
    cur = start
    pool = [st]
    stall = 0
    t0 = time.time()
    while stall <= rounds and time.time() - t0 < time_cap:
        # Same neighbour order as the substitution-based walk, so that the
        # seeded shuffle picks the same candidates.
        neighbors = []
        for st in pool:
            for i in walk:
                neighbors.append(tuple(tuple(e[:i] + (-e[i],) + e[i+1:] for e in exps) for exps in st))     # x_i -> 1/x_i
                for j in walk:
                    if i == j:
                        continue
                    neighbors.append(tuple(tuple(e[:j] + (e[j] + e[i],) + e[j+1:] for e in exps) for exps in st))   # x_i -> x_i*x_j
                    neighbors.append(tuple(tuple(e[:j] + (e[j] - e[i],) + e[j+1:] for e in exps) for exps in st))  # x_i -> x_i/x_j
        degs = [_denominator_degree(st, torus) for st in neighbors]
        m = min(degs)
        if m < cur:
            stall = 0
        else:
            stall += 1
        if m <= cur:
            cur = min(cur, m)
            cand = [st for st, dd in zip(neighbors, degs) if dd <= cur]
            rng.shuffle(cand)
            pool = cand[:concurrent]
    if cur >= start:
        return R
    verbose(f"minimize_denom_deg: denominator degree reduced from {start} to {cur}.", level=1)

    # Build h'/prod(xgens) = x^d * P'/Q', where P', Q' are the transformed
    # numerator/denominator with their monomial content removed, and
    # d = min(P) - min(Q) - torus. Put x^max(d,0) up and x^max(-d,0) down.
    num_exps, den_exps = pool[0]
    num_min, den_min = _min_exponent(num_exps), _min_exponent(den_exps)
    d = [a - b - c for a, b, c in zip(num_min, den_min, torus)]
    num_shift = [max(dk, 0) - mk for dk, mk in zip(d, num_min)]
    den_shift = [max(-dk, 0) - mk for dk, mk in zip(d, den_min)]
    num = A({tuple(ek + sk for ek, sk in zip(e, num_shift)): c
             for e, (_, c) in zip(num_exps, num_coeffs)})
    den = A({tuple(ek + sk for ek, sk in zip(e, den_shift)): c
             for e, (_, c) in zip(den_exps, den_coeffs)})
    return B(num) / B(den)


def _min_exponent(exps):
    r"""Componentwise minimum of a nonempty tuple of exponent vectors."""
    return tuple(min(col) for col in zip(*exps))


def _denominator_degree(st, torus):
    r"""
    Total degree of the reduced denominator of $h/\prod_{torus} x_i$, where
    ``st = (num_exps, den_exps)`` are the (Laurent) supports of a coprime
    numerator/denominator pair $P, Q$ of $h$. Writing $P = x^{m_P} P'$ and
    $Q = x^{m_Q} Q'$ with $P', Q'$ coprime polynomials free of monomial
    factors, the reduced denominator is $Q' x^{\max(0,\, torus - m_P + m_Q)}$.
    """
    num_exps, den_exps = st
    num_min, den_min = _min_exponent(num_exps), _min_exponent(den_exps)
    deg_q = max(sum(e) for e in den_exps) - sum(den_min)
    return deg_q + sum(max(0, c - a + b) for a, b, c in zip(num_min, den_min, torus))



def probe_reduction_order(a, f, ensure_termination, ncpus=None, profile=None, coefficient_data=None):
    r"""
    Probe r=1 and r=2 concurrently in separate single-process workers.
    Allow at most 30 seconds for r=2; r=1 has no timeout. Both start at
    the same prime, with independent bad-prime retries.
    Return ``(reduction_order, seed, degree_hint)``,
    where ``reduction_order`` is optimal choice of ``r`` and
    ``seed`` is a dict ``{prime: coeffs}`` for an operator modulo one prime.
    If ``r = 2`` is chosen, the ``ProfileData`` object ``profile`` (if given)
    is overwritten in place with the profile gathered by the ``r = 2`` probe;
    otherwise it is left unchanged. The returned degree hint schedules later
    reconstruction attempts. Only staircase entries used by the winning
    worker's successful attempt are imported into the parent's cache. Optional
    ``coefficient_data`` supplies precomputed characteristic-zero coefficient
    tables; otherwise they are prepared once here and shared with both workers.
    """
    if coefficient_data is None:
        coefficient_data = prepare_coefficient_data(a, f)
    context = multiprocessing.get_context()
    first_prime = random_prime(2**26 - 1, proof=True, lbound=2**23)

    # Separate pools let us stop r=2 without cancelling the ongoing r=1 job.
    # The OS can schedule the two single-process workers on different cores.
    with context.Pool(processes=1) as pool_r1:
        with context.Pool(processes=1) as pool_r2:
            pending_r2 = pool_r2.apply_async(
                _probe_reduction_order_worker,
                (a, f, 2, ensure_termination, first_prime, ProfileData(), 1, coefficient_data),
            )
            pending_r1 = pool_r1.apply_async(
                _probe_reduction_order_worker,
                (a, f, 1, ensure_termination, first_prime, None, 1, coefficient_data),
            )
            try:
                result_r2 = _probe_worker_result(pending_r2, timeout=30.0)
            except multiprocessing.TimeoutError:
                verbose(
                    "    Probe at reduction order 2 exceeded 30 seconds; "
                    "selecting reduction order 1.",
                    level=1,
                )
                result_r2 = None
        # Leaving the inner context terminates/joins r=2, even on timeout.
        # r=1 has been running throughout and remains alive in its own pool.
        result_r1 = _probe_worker_result(pending_r1)

    usable = [
        result
        for result in (result_r2, result_r1)
        if result is not None
    ]

    if not usable:
        # Both r=1 and r=2 were bad; start main loop at r=3, with no seed data
        # (and a disabled profile, since profiling applies only at r=2).
        return 3, {}, None

    choice = min(usable, key=lambda result: result["score"])

    seed = (
        {choice["prime"]: choice["coeffs"]}
        if choice["coeffs"] is not None
        else {}
    )

    # Adopt the r=2 profile only if r=2 is chosen; otherwise leave it disabled.
    # (ProfileData uses slots, so we copy fields explicitly rather than via __dict__.)
    if profile is not None and choice["r"] == 2:
        profile.copy_from(choice["profile"])

    # Only the winning worker's used, coefficient-independent entries survive.
    _staircase_cache.update(choice["staircases"])
    return choice["r"], seed, choice["degree_hint"]

# ~~ BOILERPLATE ~~
def _probe_reduction_order_worker(*args):
    """Transport errors whose constructors require more than Exception.args."""
    try:
        return True, _probe_reduction_order_helper(*args)
    except Exception as exc:
        return False, (type(exc), exc.args, exc.__dict__, traceback.format_exc())


def _probe_worker_result(pending, timeout=None):
    """Unpack a worker result while preserving typed errors and their fields."""
    success, result = pending.get(timeout=timeout)
    if success:
        return result
    error_type, args, attributes, worker_traceback = result
    error = error_type.__new__(error_type, *args)
    BaseException.__init__(error, *args)
    error.__dict__.update(attributes)
    error.add_note("Reduction-order probe worker traceback:\n" + worker_traceback)
    raise error


def _probe_reduction_order_helper(a, f, r, ensure_termination, first_prime=None, profile=None, ncpus=1, coefficient_data=None):
    r"""Compute one modular probe and its comparison score."""
    p = first_prime
    if coefficient_data is None:
        coefficient_data = prepare_coefficient_data(a, f)

    for attempt in range(5):
        if p is None:
            p = random_prime(2**26 - 1,proof=True,lbound=2**23,)

        try:
            # NOTE: Modifies profile in place
            staircase_keys = set()
            rho0, M, B = compute_gauss_manin_connection(
                a, f, r, p, profile=profile, staircase_keys=staircase_keys,
                coefficient_data=coefficient_data)
            
            if len(M) == 0 or rho0.nrows() == 0:
                # The reduced class vanishes: the operator is 1.
                coeffs = None
                score = (0, -1, 0, len(M), r)
            else:
                coeffs = compute_reductions_dependency(rho0, B,ensure_termination)
                degrees = [
                    -1 if c == 0 else int(c.degree())
                    for c in coeffs
                ]

                score = (
                    len(coeffs) - 1,                        # operator order
                    max(degrees),                           # maximum t-degree
                    sum(d + 1 for d in degrees if d >= 0),  # dense coefficient size
                    len(M),                                 # Stage-1 basis size
                    r,                                      # prefer r=1 on exact ties
                )

            verbose(
                f"    Probe at reduction order {r}: "
                f"operator order {score[0]}, degree {score[1]}, "
                f"basis size {len(M)}.",
                level=1,
            )

            return {"r": r, "prime": p, "coeffs": coeffs, "score": score,
                    "profile": profile,
                    "degree_hint": reconstruction_degree_hint(rho0, B),
                    "staircases": {key: _staircase_cache[key] for key in staircase_keys}}

        except Exception as exc:

            if isinstance(exc,BadPrimeError) and attempt < 4:
                p = None
                continue

            if isinstance(exc,ReductionOrderTooSmallError):
                verbose(
                    f"    Probe at reduction order {r} was unusable: "
                    f"{exc}.",
                    level=1,
                )
                return None

            raise


