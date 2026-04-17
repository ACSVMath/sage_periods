r"""
Methods and classes supporting reductions modulo the Jacobian ideal and higher order relations.
"""

# Sage imports
from sage.rings.integer import Integer
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.structure.sequence import Sequence
from sage.rings.integer_ring import ZZ
from sage.matrix.constructor import matrix
from sage.structure.element import parent


from .poly_lin_alg import echelonized_basis_poly, linear_normal_form_p

# For RhamKoszulData object.
from dataclasses import dataclass, field 

# For doing batch reductions mod a Groebner basis in Singular
# Suppress SB Warning 
from sage.libs.singular.function_factory import ff
from sage.libs.singular.option import opt_verb
opt_verb['notWarnSB'] = True

from sage.misc.verbose import verbose, set_verbose

'''
 Classes used to store (and cache) results of computations.
'''

@dataclass(slots=True)
class RhamKoszulData:
    r"""
    Container for the encoded Rham-Koszul reduction data attached to a homogeneous polynomial $f$.

    Instances of this class store the auxiliary polynomial ring with variables
    $u$ and $v$, the encoded Jacobian Gröbner basis implementing the
    Griffiths-Dwork step, the syzygy data used in the higher reductions, and
    caches for the degree-sliced spaces $U^r_q$ and $W_{\downarrow}$.
    """
    f: object
    ring: object
    xring: object | None = None
    tox: object | None = None
    fromx: object | None = None
    deg: int | Integer = 0
    df: list = field(default_factory=list)
    tsyzlm: list = field(default_factory=list)
    syz: list = field(default_factory=list)
    dim: int | Integer = 0
    jac: object | None = None
    basisWdown: dict = field(default_factory=dict)
    basis_U: dict = field(default_factory=dict)
    vars: list = field(default_factory=list)
    r: int | Integer = 1
    gens: list = field(default_factory=list)
    u: object | None = None
    v: object | None = None
    vshift: int = 0 
    # Eliminated fields variant and repmode.
    emod: object | None = None
    basis: object | None = None

    def __init__(self,f,r=1):
        r"""
        Initialise the Rham-Koszul reduction context attached to a homogeneous polynomial $f$.

        This is the setup phase for the reductions in Sections 4.3 and 7 of Lairez 2016.
        It builds the encoded polynomial ring with auxiliary variables $u$, $v$,
        prepares the Jacobian ideal implementing the Griffiths-Dwork step, computes the
        relevant syzygy data, and initialises the caches that will later store the
        degree-sliced spaces $U^r_q$ and $W_{\downarrow}$.

        INPUT:

        * ``f`` -- A homogeneous multivariate polynomial.
        * ``r`` -- The reduction order to be used later for ``[]_r``.

        OUTPUT:

        * ``U`` -- A record ``U`` of type ``RhamKoszulDataformat`` containing all data
          needed by ``hom_reduce``, ``gauss_manin_helper``, and the other reduction routines.
        """
    
        self.f=f
        self.r=r
        self.dim=parent(f).ngens()
        self.ring=parent(f)
        self.basis_U={}
        self.basisWdown={}
        self.deg=f.degree()
        
        n = self.dim
        names = ["u", "v"] + [f"x{i}" for i in range(n)]
        R = PolynomialRing(self.ring.base_ring(), n + 2, names=names, order="degrevlex")
    
        # Default encoding: variables are ordered u, v, x0, ..., x_{n-1}.  The maps
        # tox and fromx move between the original polynomial ring and the extended one.
        u = R.gen(0)
        v = R.gen(1)
        tox = self.ring.hom([R.gen(i+2) for i in range(n)],R)  # to extended ring
        fromx = R.hom([self.ring.one(), self.ring.zero()] + [self.ring.gen(i) for i in range(n)],self.ring)  # from extended ring
        self.vshift = 2
    
        self.u = u
        self.v = v
        self.xring = R
        self.tox = tox
        self.fromx = fromx
    
        # Store df and the x-variables in the encoded ring
        self.df = [tox(f.derivative(self.ring.gen(i))) for i in range(n)]
        self.vars = [tox(self.ring.gen(i)) for i in range(n)]
    
        # The powers u^(n-i) v^i encode the standard basis used to represent the form
        # module inside the polynomial ring A[u,v]. See the discussion in Section 7.1 of Lairez
        # on the construction of the polynomial module M'.
        self.gens = [u ** (n - i) * v ** i for i in range(n + 1)]
        xrels = [u ** i * v ** (n + 1 - i) for i in range(n + 2)]
    
        # Build the ideal implementing the Griffiths-Dwork reduction step in encoded
        # form: the relations df_i * omega - xi_i, together with the monomial relations
        # cutting back to the encoded module.

        # Allow passing to Singular and batching based on whether or not r > 1.
        # (Note that in the smooth case, it is a bad idea to pass to singular, 
        # because we add the overhead of converting to and from Singular with 
        # no gains from batching.)
        if self.r > 1:
            self.jac = ff.groebner(self.xring.ideal([self.df[i] * self.gens[0] - self.gens[i + 1] for i in range(n)] + xrels))
        else:
            self.jac = self.xring.ideal([self.df[i] * self.gens[0] - self.gens[i + 1] for i in range(n)] + xrels).groebner_basis()
        # We only keep the Groebner basis. The original generators can be found in the above way, easily.
    
        # Compute Groebner basis for Syz; this suffices, via Lemma 32 in Lairez 2016
        self.syz = [p for p in self.jac if p.degree(u) < n]

        # Generate Groebner basis for lms of Trivial Syzygies
        triv_syz = [self.df[i] * self.gens[j] - self.df[j] * self.gens[i] for i in range(n) for j in range (n)] + xrels
        self.tsyzlm = list(self.xring.ideal([p.lm() for p in triv_syz]).groebner_basis()) if triv_syz else []

@dataclass(slots=True)
class RKGaussManinData:
    r"""
    Container for the finite-dimensional data returned by ``gauss_manin_helper``.

    Its fields store the monomial basis closed under the reduced connection action,
    the corresponding exponent vectors, the projection of the starting classes onto
    that basis, and the matrix of multiplication by ``-der`` followed by reduction.
    """
    basis: list = field(default_factory=list)
    ebasis: list = field(default_factory=list)
    proj: object | None = None
    gm: object | None = None


'''
 These helpers build degree slices of syzygy spaces.  They work entirely 
 at the level of monomials and leading monomials, so that the costly 
 linear-algebra part is replaced by staircase computations and Groebner reductions.
'''

def _diff_stairs(a, B, vars, deg):
    r"""
    Compute the degree ``deg`` staircase descendants of the monomial ``a`` modulo ``B``.

    This is a degree-by-degree staircase walk. Starting from ``a``, it multiplies by
    the variables in ``vars`` until total degree ``deg`` is reached, reducing modulo the
    forbidden monomials ``B`` at each stage. The result is a monomial basis for the
    relevant degree slice outside the leading-monomial ideal generated by ``B``.

    INPUT:

    * ``a`` -- A monomial, represented as an element of a Sage multivariate polynomial ring.
    * ``B`` -- A list of polynomials or monomials in the same ring, forming a Groebner basis, 
      used through their leading monomials as the forbidden staircase boundary.
    * ``vars`` -- A list of generators of ``parent(a)`` used to build staircase descendants.
    * ``deg`` -- A non-negative integer specifying the target total degree.

    OUTPUT:

    * ``descendants`` -- A Python ``set`` of monomials in ``parent(a)`` giving the
      admissible degree ``deg`` descendants of ``a`` modulo ``B``.
    """
    if deg < a.degree():
        # No descendant can have smaller total degree than the starting monomial.
        return set()
    if deg == a.degree():
        # At the target degree, keep the normal form if it does not vanish.
        nf = a.reduce(B)
        return {nf} - {a.parent().zero()}

    # Recurse one degree lower, then extend by one variable and reduce again.
    recurse = _diff_stairs(a, B, vars, deg - 1)
    ret = {v * e for v in vars for e in recurse}
    ret = set([r.reduce(B) for r in ret]) - {a.parent().zero()}
    return ret


def basis_syzygies(U, deg):
    r"""
    Return a degree slice of non-trivial syzygy space used in Algorithm 3 of Lairez 2016.

    In the notation of Lairez, the faster recursion uses a complement $A_q$ of the
    trivial syzygies $S'_q$ inside the full syzygy space $S_q$. This routine extracts
    one homogeneous slice of that complement from the Gröbner basis data prepared
    in the init method of ``RhamKoszulData``.

    INPUT:

    * ``U`` -- The reduction context record produced by initialization of
      ``RhamKoszulData`` class.
    * ``deg`` -- An integer specifying the polynomial degree slice to extract.

    OUTPUT:

    * ``syzygies`` -- A sequence of multivariate polynomials representing a basis of
      the requested homogeneous slice of non-trivial syzygies.
    """

    # Build a basis of a degree slice generated by A modulo the forbidden monomials B.

    # For each generator a in A, enumerate the admissible staircase monomials above its leading monomial
    # then lift them back to actual multiples of a. This is used to extract degree q 
    # pieces of the non-trivial syzygy space A_q appearing in the faster construction of X^r_q.
    ret = []
    BB = list(U.tsyzlm)

    for a in U.syz:
        lm = a.lm()

        # Compute the admissible degree-deg monomials divisible by lm but not by the
        # leading monomials already seen.
        D = _diff_stairs(lm, BB, U.vars, deg+1)
        ret.extend([(d // lm) * a for d in D])

        # Augment the forbidden set so that the output is echelon-like and contains
        # no duplicate leading monomial.
        BB.append(lm)
    return ret


def exterior_derivative(U, a):
    r"""
    Compute the encoded exterior differential $d$ on the auxiliary representation.

    The polynomial $a$ is interpreted in the $u, v$-encoding of forms, as a part of the
    ``RhamKoszulData`` structure $U$. Splitting $a$ by powers of $v$ recovers the coefficients
    of the form components, and the returned polynomial is the encoded image under $d$.

    This is an internal function for sage_periods, and is not 
    meant to be called by the user.

    INPUT:

    * ``U`` -- A ``RhamKoszulData`` structure.
    * ``a`` -- A polynomial, interpreted in the ``u``, ``v`` encoding of forms.

    OUTPUT:

    * ``da`` -- A polynomial representing the top-form $d(a)$ in our encoded representation.

    ASSUMPTIONS:

    * The polynomial ``a`` represents an n-form in our encoded ring. This will not work if other-grade forms are present!

    EXAMPLES:

    A nontrivial example.

            sage: R = PolynomialRing(QQ, ["x0", "x1"])
            sage: x0, x1 = R.gens()
            sage: f = x0**2 + x1**2
            sage: U = RhamKoszulData(f)
            sage: X0, X1 = U.vars
            sage: a = 3*X0*U.u*U.v + 2*X1*U.v**2 # This corresponds to the form \alpha = 3x0xi[0]  +2x1 xi[1]
            sage: exterior_derivative(U, a)
            5*u^2
            
    A trivial example.

            sage: R = PolynomialRing(QQ, ["x0", "x1"])
            sage: x0, x1 = R.gens()
            sage: f = x0**2 + x1**2
            sage: U = RhamKoszulData(f)
            sage: exterior_derivative(U, U.xring.zero())
            0 

    This function is meant to be applied to $n$-forms. If you apply it
    to, for instance, an $(n+1)$-form, it will not work correctly.

            sage: eta = 3*U.u**2 
            sage: # In form-land the output should be 0, but this will return eta itself.
            sage: exterior_derivative(U,eta)
            3*u^2      

    """

    c = [a.coefficient({U.v:i}) for i in range(a.degree(U.v)+1)]
    if len(c) > 0:
        sh = U.vshift
        return c[0] + sum(
            c[i].derivative(parent(c[i]).gen(i - 1 + sh)) * U.u ** i
            for i in range(1, len(c))
        )
    return U.xring.zero()


def elementary_reduction_step(U, L0):
    r"""
    Perform one batched Griffiths-Dwork reduction step on the encoded family ``L0``.

    Each element of ``L0`` is first reduced modulo the encoded Jacobian Gröbner basis
    stored in ``U.jac``. The remaining $\beta$-part is then sent through
    ``exterior_derivative``, producing the next batch that appears in the recursive
    construction of the higher reductions.

    This is an internal function for sage_periods, and is not 
    meant to be called by the user.

    INPUT:

    * ``U`` -- A ``RhamKoszulData`` object containing the encoded Jacobian data.
    * ``L0`` -- A Python list of encoded polynomials in ``U.xring`` to be reduced in batch.

    OUTPUT:

    * ``L`` -- A Python list of encoded polynomials obtained after one batched
      Griffiths-Dwork reduction step.

    EXAMPLES:

    A nontrivial example.

            sage: R = PolynomialRing(QQ, ["x0", "x1"])
            sage: x0, x1 = R.gens()
            sage: f = x0**2 + x1**2
            sage: U = RhamKoszulData(f)
            sage: X0, X1 = U.vars
            sage: L0 = [X0*U.v, X1*U.v**2]
            sage: U.jac
            [u^3, u^2*v, u*v^2, v^3, u^2*x0 - 1/2*u*v, v^2*x0 - u*v*x1, u^2*x1 - 1/2*v^2]
            sage: L0[0].reduce(U.jac)
            v*x0
            sage: L0[1].reduce(U.jac)
            v^2*x1
            sage: # These two encoded polynomials are already in normal form
            sage: # modulo the Grobner basis of Jac, and exterior_derivative sends
            sage: # x0*v to u and x1*v^2 to u^2.
            sage: elementary_reduction_step(U, L0)
            [u, u^2]
            
    A trivial example.

            sage: R = PolynomialRing(QQ, ["x0", "x1"])
            sage: x0, x1 = R.gens()
            sage: f = x0**2 + x1**2
            sage: U = RhamKoszulData(f)
            sage: elementary_reduction_step(U, [])
            []

    """
    if len(L0) == 0:
        return L0

    if U.r >1: # Should we batch in Singular, or just do 1-by-1 in sage?
        L = ff.reduce(U.xring.ideal(L0),U.jac)
    else:
        L = [p.reduce(U.jac) for p in L0]
    return [exterior_derivative(U, p) for p in L]


'''
  The following routines implement the fast recursive higher reductions.  The
  code works degree by degree and stores the spaces
  in a split form adapted to batched echelonisation.
'''


def basis_U(U, r, q):
    r"""
    Compute the degree $q$ slice $U^r_q$ used by the fast higher reduction $[]_r$.

    This is a code-level refinement of the spaces $X^r_q$ from Section 4.3 of Lairez 2016.
    Instead of storing the whole $X^r_q$ at once, the implementation keeps:

    - ``basis_U(r, q)``: the exact degree ``q`` pivot part, and
    - ``basisWdown(r, q)``: the strictly lower-degree tail (called $W_{\downarrow}$ here).

    Together these two lists encode the same information needed from $X^r_q$ for the
    recursive reduction.

    This is an internal function for sage_periods, and is not meant to be called by the user.

    INPUT:

    * ``U`` -- An ``RhamKoszulData`` object, holding our Rham-Koszul complex / reduction data.
    * ``r`` -- A non-negative integer specifying the order of the higher reduction.
    * ``q`` -- An integer specifying the degree of the slice taken.

    OUTPUT:

    * ``None`` -- No direct return value. The objects ``U.basis_U[(r, q)]`` and, when relevant,
      ``U.basisWdown[(r, q)]`` are filled in place.

    EXAMPLES:

    A nontrivial example.

            sage: R = PolynomialRing(QQ, ["x0", "x1"])
            sage: x0, x1 = R.gens()
            sage: f = x0**2 + x1**2
            sage: U = RhamKoszulData(f)
            sage: # For this f, Jac(f) is 0-dimensional, so V(f) is smooth.
            sage: # Thus all syzygies are trivial syzygies, so they will not be apart of this space.
            sage: basis_U(U, 1, 4)
            sage: U.basis_U[(1, 4)]
            []
            sage: U.basisWdown[(1, 4)]
            [0]
            

    A trivial example.

            sage: R = PolynomialRing(QQ, ["x0", "x1"])
            sage: x0, x1 = R.gens()
            sage: f = x0^2 + x1^2
            sage: U = RhamKoszulData(f)
            sage: basis_U(U, 0, 4)
            []

    !!! note

        Technically [0] is not a basis, but this doesn't slow down futher computations.

    """
    
    if r < 0:
        raise ValueError("r must be non-negative")
    if (r, q) in U.basis_U:
        # Cached value already available.
        return None
    if r == 0:
        # There is no higher reduction space at level 0.
        U.basis_U[(r, q)] = []
        return None
    if r == 1:
        # Base case: X^1_q = d A_{q-1}.  The non-trivial syzygies supply the A-part,
        # and cast_syzygies applies d in the chosen representation mode.
        syz = basis_syzygies(U, q - U.deg)
        U.basisWdown[(r, q)] = [exterior_derivative(U, p) for p in syz] # cast_syzygies, in the r=1 case.
        U.basis_U[(r, q)] = []
        return None
    if r > 1:
        # Recursive step: use the already-computed data at order r-1, both in degree
        # q and in degree q + deg(f).  This mirrors the recursion
        # X^{r}_q = dA_{q-1} + red^GD_q( X^{r-1}_{q+deg(f)} cap F_q ).
        basis_U(U, r - 1, q)
        basis_U(U, r - 1, q + U.deg)

        # Push the stored lower-degree tail down by one Griffiths-Dwork step.
        rels = elementary_reduction_step(U, U.basisWdown[(r - 1, q + U.deg)])

        # Echelonize together the previous degree-q slice and the newly produced relations.  
        new = echelonized_basis_poly(U.basis_U[(r-1,q)]+rels, U.xring, return_matrix=False)

        # Split the resulting echelon basis into the exact-degree-q pivot part and the strict lower-degree tail.
        U.basis_U[(r, q)], U.basisWdown[(r, q)] = [p for p in new if p.degree() == q], [p for p in new if p.degree() < q]
        return None


def _hom_reduce_helper(U, L, r):
    r"""
    Recursive worker for ``hom_reduce``.

    Starting from the current top degree ``q``, the function repeatedly performs one batched GD step,
    then eliminates the degree-``q`` component against the precomputed space $U^r_q$.

    INPUT:

    * ``U`` -- A ``RhamKoszulData`` object carrying the cached reduction data.
    * ``L`` -- A mutable Python list of encoded polynomials in ``U.xring``. It is modified in place.
    * ``r`` -- A positive integer specifying the order of the higher reduction.

    OUTPUT:

    * ``None`` -- No direct return value. The list ``L`` is rewritten in place by its reduced encoded representatives.
    """
    if len(L) == 0:
        return None

    q = ZZ(max(p.degree() for p in L)) # compute_weight
    while q > 0:
        # First apply one Griffiths-Dwork reduction step to the whole batch.
        newL = elementary_reduction_step(U, list(L))
        L[:] = newL

        if r > 1:
            # Then quotient out by the precomputed degree-q part U^r_q.
            basis_U(U, r, q)
            L[:] = linear_normal_form_p(list(L), U.basis_U[(r, q)], False)
            '''
                TODO: It *might* be faster to write a version of linear_normal_form
                that works in some VectorSpace / CombinatorialFreeModule,
                rather than working with polynomials naively.
            '''

        # Each step lowers the relevant degree by deg(f).
        q -= U.deg
    return None


def hom_reduce(U, L, r):
    r"""
    Reduce a list of homogeneous numerators modulo the higher reduction $[\cdot]_r$.

    This is the user-facing wrapper around the internal encoded recursion. It
    converts ordinary polynomials to the internal top-form representation, applies
    the recursive reduction based on the precomputed $U^r_q$ spaces, and decodes the
    result back to the original polynomial ring.

    This is an internal function for sage_periods, and is not 
    meant to be called by the user.

    INPUT:

    * ``U`` -- A ``RhamKoszulData`` object, passed "by reference" for cache reuse.
    * ``L`` -- A sequence of polynomials in ``U.ring`` of the appropriate homogeneous degree.
    * ``r`` -- A non-negative integer representing the order of the reduction.

    OUTPUT:

    * ``None`` -- No direct return value. The list ``L`` is replaced in place by its
      ``[]_r``-reduced representatives.

    EXAMPLES:

    A nontrivial example.

            sage: R = PolynomialRing(QQ, ["x0", "x1"])
            sage: x0, x1 = R.gens()
            sage: f = x0**2 + x1**2
            sage: U = RhamKoszulData(f)
            sage: L = [R(1), x0, x1]
            sage: # The terms x0 and x1 lie in the Jacobian ideal generated by 2*x0 and 2*x1 
            sage: # and therefore reduce to 0, while the constant class 1 is unchanged.
            sage: hom_reduce(U, L, 1)
            [1, 0, 0]
            
    A trivial example.

            sage: R = PolynomialRing(QQ, ["x0", "x1"])
            sage: x0, x1 = R.gens()
            sage: f = x0**2 + x1**2
            sage: U = RhamKoszulData(f)
            sage: L = []
            sage: hom_reduce(U, L, 1)
            []
    """

    if r < 1:
        raise ValueError("r must be positive")
    if len(L) == 0:
        return None
    if Sequence(list(L)).universe() != U.ring:
        raise TypeError("Bad argument type")

    # Switch to the encoded model, reduce there, then decode back.
    enc = [U.tox(p) * U.gens[0] for p in list(L)] 
    _hom_reduce_helper(U, enc, r)
    L[:] = [U.fromx(p) for p in enc] 

    return None


def gauss_manin_helper(U, der, L, ret=None):
    r"""
    Compute the finite-dimensional linear action induced by multiplication by ``-der``.

    One works on the smallest monomially generated subspace stable under
    the reduction of the connection operator. Starting from the classes represented
    by $L$, this routine constructs that monomial basis frontier-by-frontier, reduces
    all images in batch, and returns both the action matrix and the projection of $L$
    onto the resulting basis.

    This is an internal function for sage_periods, and is not 
    meant to be called by the user.

    INPUT:

    * ``U`` -- An ``RhamKoszulData`` object, passed "by reference" for cache reuse.
    * ``der`` -- A polynomial in ``U.ring``, typically the evaluated $f^\delta$.
    * ``L`` -- A sequence of polynomials whose reduced classes generate the starting space.
    * ``ret`` -- An optional output record variable passed by reference.

    OUTPUT:

    * ``ret`` -- An ``RKGaussManinData`` record whose fields are filled with:
        * ``basis`` -- The monomial basis found by the closure process.
        * ``ebasis`` -- The exponent vectors of ``basis``.
        * ``proj`` -- The matrix expressing the reduced input classes in that basis.
        * ``gm`` -- The matrix of the map $P \mapsto [-\mathrm{der} P]_r$ in that basis.

    EXAMPLES:

    A nontrivial example.

            sage: R = PolynomialRing(QQ, ["x0", "x1"])
            sage: x0, x1 = R.gens()
            sage: f = x0**2 + x1**2
            sage: U = RhamKoszulData(f)
            sage: 
            sage: # Expected outputs:
            sage: #   basis  = [1]
            sage: #   ebasis = (((0, 0),),)
            sage: #   proj   = [1]
            sage: #   gm     = [-1]
            sage: # Justification: the reduced class of 1 spans a one-dimensional space, and
            sage: # multiplication by -1 acts on it by the scalar -1.
            sage: ret = gauss_manin_helper(U, R(1), [R(1)])
            sage: ret.basis
            [1]
            sage: ret.ebasis
            (((0, 0),),)
            sage: ret.proj
            [1]
            sage: ret.gm
            [-1] 

    A trivial example.

            sage: R = PolynomialRing(QQ, ["x0", "x1"])
            sage: x0, x1 = R.gens()
            sage: f = x0**2 + x1**2
            sage: U = RhamKoszulData(f)
            sage: ret = gauss_manin_helper(U, R(1), [])
            sage: ret.basis
            []
            sage: ret.ebasis
            ()
            sage: ret.proj
            []
            sage: ret.gm
            []
            
    """

    basis = []
    gm = []

    # Reduce the initial generators and extract the monomials supporting them.
    proj = list(L)
    hom_reduce(U, proj, U.r)
    verbose("Initial proj after hom_reduce:",level=1)
    verbose(proj,level=1)
    xmons = list(set().union(*(set(p.monomials()) for p in proj))) if proj else []

    while len(xmons) > 0:
        # Add the current frontier of monomials to the basis.
        xmons = sorted(xmons)
        basis.extend(xmons)
        verbose("xmons:",level=1)
        verbose(xmons,level=1)
        verbose("basis:",level=1)
        verbose(basis,level=1)

        # Compute the images of the current frontier under the connection term and
        # reduce them modulo []_r, all in batch.
        nf = [-der * p for p in xmons]
        verbose("nf: ",level=1)
        verbose(nf,level=1)
        hom_reduce(U, nf, U.r)
        verbose("hom_reduced nf:",level=1)
        verbose(nf,level=1)
        gm.extend(nf)

        # If the reduced images exceed the expected degree bound, then r was too
        # small for the higher reduction to have stabilised.
        if max(p.degree() for p in gm) > (U.dim - 1) * U.deg:
            raise RuntimeError("INCREASE_R")

        # Continue with the newly discovered monomials, exactly as in the monomial
        # closure construction of M in Algorithm 4 of Lairez 2016.
        newmons = set().union(*(set(p.monomials()) for p in nf))
        xmons = list(newmons.difference(set(basis)))
        verbose("xmons after subtracting against basis:",level=1)
        verbose(xmons,level=1)

    ret = RKGaussManinData() if ret is None else ret
    ret.basis = basis

    # Convert the reduced images to the matrix of the action in the monomial basis.
    ret.gm = matrix(
        U.ring.base_ring(),
        len(basis),
        len(basis),
        [[gm[j].monomial_coefficient(basis[i]) for j in range(len(basis))] for i in range(len(basis))],
    )

    # Likewise, express the starting generators in that same basis.
    ret.proj = matrix(
        U.ring.base_ring(),
        len(basis),
        len(L),
        [[proj[j].monomial_coefficient(basis[i]) for j in range(len(L))] for i in range(len(basis))],
    )

    ret.ebasis = tuple(tuple(m.exponents(False)) for m in basis)
    return ret

__all__ = [
    "basis_syzygies",
    "gauss_manin_helper",
    "elementary_reduction_step",
    "exterior_derivative",
    "hom_reduce",
    "RKGaussManinData",
    "RhamKoszulData",
]