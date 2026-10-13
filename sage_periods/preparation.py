r"""Functions used in preparing rational period intetgral for modular pipeline."""

from sage.all import *

from .pipeline import compute_gauss_manin_connection, compute_reductions_dependency
from .errors import *
from .profiling import ProfileData

from sage.misc.verbose import verbose
import multiprocessing
import pickle
import tempfile
import time
import traceback

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

# TODO
def minimize_denominator_degree(R, xgens, rounds=3, concurrent=10, time_cap=30.0):
    r"""
    Torus (monomial) change of variables lowering the denominator degree.

    Port of the ``minimize_degree`` heuristic of Lairez's MAGMA package: a
    random walk over the unimodular substitutions $x_i \mapsto 1/x_i$ and
    $x_i \mapsto x_i x_j^{\pm 1}$, applied to the $\dd x/x$-normalized
    integrand $h = R\,\prod_i x_i$, keeping a pool of the best candidates
    and stopping after ``rounds`` rounds without improvement (or after
    ``time_cap`` seconds).  Returns the transformed integrand
    $h'/\prod_i x_i$, or ``R`` itself if no improvement was found.

    This is an internal function for sage_periods, and is not meant to be
    called by the user.
    """
    import random as _random
    import time as _time
    rng = _random.Random(20260803)
    prodx = prod(xgens)
    h = R * prodx

    # Restrict the walk to the variables that actually occur in h: the
    # integrand may be independent of some torus coordinates (for instance
    # the diagonal integrand of 1/(1-x*y) collapses to 1/((1-t)*y), whose
    # h = R*prod(x) is the constant 1/(1-t)). A monomial substitution in
    # an absent variable fixes h, so we needn't include them.
    Agens = {str(g): g for g in h.numerator().parent().gens()}
    walkgens = [v for v in xgens
                if h.numerator().degree(Agens[str(v)]) > 0
                or h.denominator().degree(Agens[str(v)]) > 0]
    if not walkgens:
        return R

    start = (h / prodx).denominator().degree()
    cur = start
    pool = [h]
    stall = 0
    t0 = _time.time()
    while stall <= rounds and _time.time() - t0 < time_cap:
        neighbors = []
        for hh in pool:
            for i, xi in enumerate(walkgens):
                neighbors.append(hh.subs({xi: 1 / xi}))
                for j, xj in enumerate(walkgens):
                    if i == j:
                        continue
                    neighbors.append(hh.subs({xi: xi * xj}))
                    neighbors.append(hh.subs({xi: xi / xj}))
        degs = [(hh / prodx).denominator().degree() for hh in neighbors]
        m = min(degs)
        if m < cur:
            stall = 0
        else:
            stall += 1
        if m <= cur:
            cur = min(cur, m)
            cand = [hh for hh, dd in zip(neighbors, degs) if dd <= cur]
            rng.shuffle(cand)
            pool = cand[:concurrent]
    if cur >= start:
        return R
    verbose(f"minimize_denom_deg: denominator degree reduced from {start} to {cur}.", level=1)
    return pool[0] / prodx

def _probe_reduction_order_helper(a, f, r, ensure_termination, max_basis=None):
    """Return (modular result, basis size, profile) for one reduction order."""
    profile = ProfileData() if r == 2 else None
    for attempt in range(5):
        p = random_prime(2**26 - 1, proof=True, lbound=2**23)
        try:
            rho0, M, B = compute_gauss_manin_connection(
                a, f, r, p, profile=profile, max_basis=max_basis)
            if len(M) == 0 or rho0.nrows() == 0:
                result = (0, p, None)
            else:
                coeffs = compute_reductions_dependency(rho0, B, ensure_termination)
                result = (len(coeffs) - 1, p, coeffs)
            return result, len(M), profile
        except (ReductionOrderTooSmallError, ProbeBasisCapExceededError):
            return (None, None, None), None, None
        except BadPrimeError:
            if profile is not None:
                profile.reset()
            if attempt == 4:
                raise


def _probe_worker(a, f, ensure_termination, result_path):
    """Write one complete result; no pipe/queue can block worker shutdown."""
    try:
        payload = (True, _probe_reduction_order_helper(a, f, 2, ensure_termination))
    except Exception as exc:
        # Our typed errors require more constructor arguments than exc.args.
        # Transport their state explicitly rather than pickling the exception.
        payload = (False, (type(exc), exc.args, exc.__dict__, traceback.format_exc()))
    with open(result_path, "wb") as stream:
        pickle.dump(payload, stream, protocol=pickle.HIGHEST_PROTOCOL)


def _timed_probe_reduction_order(a, f, ensure_termination, timeout=30.0):
    """Run r=2 under a single deadline, including retries and result writing.

    Return None on timeout. Cleanup may add at most two half-second joins;
    unexpected worker failures are errors, not timeouts or bad primes.
    """
    # Fork preserves the already-initialized Sage environment on Unix. Other
    # platforms use their default start method and must support Sage pickling.
    # !! NOTE: Add warning about this to compatibility in documentation!

    # ~~ BOILERPLATE ~~
    method = "fork" if "fork" in multiprocessing.get_all_start_methods() else None
    context = multiprocessing.get_context(method)
    with tempfile.TemporaryDirectory(prefix="sage-periods-probe-") as directory:
        result_path = directory + "/result.pkl"
        worker = context.Process(target=_probe_worker, args=(a, f, ensure_termination, result_path))
        deadline = time.monotonic() + timeout
        try:
            worker.start()
            worker.join(max(0.0, deadline - time.monotonic()))
            if worker.is_alive():
                return None
            if worker.exitcode != 0:
                raise RuntimeError(f"Reduction-order probe worker exited with code {worker.exitcode}.")
            with open(result_path, "rb") as stream:
                success, payload = pickle.load(stream)
            if success:
                return payload
            error_type, args, attributes, worker_traceback = payload
            error = error_type.__new__(error_type, *args)
            BaseException.__init__(error, *args)
            error.__dict__.update(attributes)
            error.add_note("Reduction-order probe worker traceback:\n" + worker_traceback)
            raise error
        finally:
            if worker.is_alive():
                worker.terminate()
                worker.join(0.5)
                if worker.is_alive():
                    worker.kill()
                    worker.join(0.5)
                if worker.is_alive():
                    raise RuntimeError("Could not stop the reduction-order probe worker.")
            worker.close()


def probe_reduction_order(a, f, ensure_termination=False, max_probe_order=2,
                          ncpus=None, probe_timeout=30.0):
    r"""
    Probe r=2 in a timed worker, then r=1 with a basis cap, then higher orders
    up to ``max_probe_order``. Return ``(r, seed, profile)``.

    Each successful probe stores ``(operator_order, prime, coeffs)``. Select
    the smallest operator order, breaking ties in favour of the smaller r.
    An operator of order <= 1 triggers an early exit, as a selection heuristic.
    If no probe is usable, start the main loop at ``max_probe_order + 1``.

    ``probe_timeout`` is one deadline for the entire r=2 probe, including
    Stage 1, Stage 2 and all bad-prime retries. On timeout its partial data
    is discarded and r=1 is tried without a cap. On success, r=1 is capped
    at ``3*max(1, size_at_r2) + 8``. Other orders and the main computation
    have no wall-clock timeout. Worker cleanup can briefly extend the wait.

    ``seed`` is ``{prime: coeffs}``, or empty for a vanishing reduced class
    or when no probe succeeds. ``profile`` is a ``ProfileData`` object only
    when r=2 is selected; otherwise it is ``None``. Bad primes are retried
    up to five attempts, clearing any profile from the discarded attempt.
    Unexpected errors propagate, including errors raised inside the worker.

    ``ensure_termination`` selects the exact Stage-2 rank calculation.
    ``ncpus`` is retained for API compatibility; probes run sequentially.
    """

    results = {}
    profiles = {}
    basis_sizes = {}
    probe_orders = (list(range(min(2, max_probe_order), 0, -1))
                    + list(range(3, max_probe_order + 1)))
    
    for r in probe_orders:
        cap = None
        if r == 1 and 2 in basis_sizes:
            cap = 3 * max(1, basis_sizes[2]) + 8 # TODO: Find better heuristic than this?
        if r == 2:
            outcome = _timed_probe_reduction_order(a, f, ensure_termination, probe_timeout)
            if outcome is None:
                verbose(f"    Probe at reduction order 2 exceeded {probe_timeout:g} seconds; "
                        "trying reduction order 1.", level=1)
                results[r] = (None, None, None)
                continue
        else:
            outcome = _probe_reduction_order_helper(a, f, r, ensure_termination, cap)
        results[r], basis_size, profile = outcome
        if basis_size is not None:
            basis_sizes[r] = basis_size
        if profile is not None:
            profiles[r] = profile

        order = results[r][0]
        verbose(f"    Probe at reduction order {r}: modular operator order {order}.", level=1)
        if order is not None and order <= 1:
            break

    usable = [(result[0], r) for r, result in results.items() if result[0] is not None]
    if not usable:
        return max_probe_order + 1, {}, None
    _, choice = min(usable)
    _, p, coeffs = results[choice]
    seed = {p: coeffs} if coeffs is not None else {}
    return choice, seed, profiles.get(choice)
