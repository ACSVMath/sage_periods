r"""
Linearized reduction closure over prime fields.

This module reimplements the inner loop of the plain (non-certified)
higher reduction -- the Griffiths-Dwork step, the recursive spaces
$U^r_q$ / $W_\downarrow$, and the monomial closure of
``gauss_manin_helper`` -- as linear algebra over an indexed monomial
basis, instead of polynomial arithmetic term by term.

Everything the original code does to a polynomial is a linear operation
on its coefficient vector:

* the Singular normal form modulo the fixed Groebner basis ``U.jac`` is
  a linear projection (the reduced normal form modulo a fixed Groebner
  basis is unique, hence linear in the input), so it is computed once
  per *monomial* per evaluation point and applied to whole batches by
  substitution;
* the encoded exterior differential sends the monomial
  $u^a v^i x^E$ (for $i \ge 1$) to $E_{i-1}\, u^{a+i} x^{E - e_{i-1}}$
  and fixes $v$-free monomials, so it is a precomputed scaled
  column-scatter;
* echelonization and elimination against an echelon basis act on
  coefficient rows directly.

A batch element is stored as a dict ``{column: coefficient}`` with
integer coefficients in ``[0, p)``; columns index encoded-ring
monomials in order of first appearance, with their degrevlex sort keys
cached.  Dense echelonization goes through the same
``Matrix(GF(p), ...).echelonize()`` kernel as before, with the active
columns sorted by the same degrevlex key, so pivot choices -- and hence
every output -- coincide exactly with the polynomial implementation.

The engine is used automatically by ``gauss_manin_helper`` when the
coefficient field is a prime finite field (i.e., in the modular pipeline);
the polynomial implementation remains in place for other base rings and 
for the certificate variant.
"""

import numpy as np
import scipy.sparse as _sp

from sage.matrix.constructor import Matrix, matrix
from sage.libs.singular.function_factory import ff
from sage.rings.polynomial.polydict import ETuple


class LinearEngine:
    r"""
    Per-evaluation-point linearization context.

    Holds the growing monomial/column tables of the encoded ring of a
    ``RhamKoszulData`` instance, the per-monomial normal-form memo, the
    exterior-derivative scatter tables, and the engine-side caches for
    the degree slices $U^r_q$ and $W_\downarrow$ (stored as coefficient
    rows rather than polynomials).

    This is an internal class for sage_periods, and is not meant to be
    used directly.
    """

    def __init__(self, U):
        self.U = U
        self.p = int(U.ring.base_ring().characteristic())
        self.vshift = U.vshift
        # Column tables for the encoded ring: exponent tuple <-> column.
        self.exps = []       # column -> exponent tuple
        self.col_of = {}     # exponent tuple -> column
        self.deg = []        # column -> total degree
        self.skey = []       # column -> degrevlex sort key
        self.ext_tgt = []    # column -> exterior derivative target column
                             #   (self for v-free columns, -1 for kernel)
        self.ext_mul = []    # column -> exterior derivative multiplier
        self.nf = []         # column -> None (unknown) | True (normal
                             #   monomial) | row dict (its normal form)
        self._skey_fn = U.xring.term_order().sortkey
        jac_lms = [tuple(g.lm().exponents()[0]) for g in U.jac]
        # Dense integer matrix of the Groebner leading exponents, for the
        # vectorized normality (staircase membership) test.
        self._jaclm_arr = (np.array(jac_lms, dtype=np.int64)
                           if jac_lms else None)
        # Column tables for the decoded ring U.ring (x-variables only).
        self._xskey_fn = U.ring.term_order().sortkey
        self._xskey_cache = {}
        # Engine-side caches mirroring U.basis_U / U.basisWdown.
        # basis_U values are lists of (row, pivot_column); basisWdown
        # values are lists of rows.
        self.basis_U = {}
        self.basisWdown = {}
        self.wdown_liftlms = {}
        # Per-column composed Griffiths-Dwork images d(NF(c)), stored as
        # immutable (target_columns, values) tuples.  These are fixed once
        # computed because both the canonical normal form and exterior-
        # derivative scatter are fixed for the lifetime of the engine.
        self.step_arr = {}
        # Coefficient rows of the non-trivial syzygy generators, used to
        # build staircase lifts by exponent shifts.
        self.syz_terms = [
            [(tuple(e), int(c)) for e, c in s.dict().items()] for s in U.syz
        ]

    # ------------------------------------------------------------------
    # column bookkeeping

    def col(self, e):
        r"""
        Return the column of the encoded exponent tuple ``e``, creating it
        (with degree, sort key, exterior-derivative image and an unknown
        normal form) on first sight.
        """
        c = self.col_of.get(e)
        if c is not None:
            return c
        # If not cache, build.
        c = len(self.exps)
        self.col_of[e] = c
        self.exps.append(e)
        self.deg.append(sum(e))
        self.skey.append(self._skey_fn(ETuple(list(e))))
        self.nf.append(None)
        # Exterior derivative action on this monomial.
        self.ext_tgt.append(0)
        self.ext_mul.append(0)
        i = e[1]
        if i == 0:
            self.ext_tgt[c] = c
            self.ext_mul[c] = 1
        else:
            xidx = i - 1 + self.vshift
            k = e[xidx]
            if k == 0:
                self.ext_tgt[c] = -1
            else:
                enew = list(e)
                enew[1] = 0
                enew[0] = e[0] + i
                enew[xidx] = k - 1
                tc = self.col(tuple(enew))  # depth-one recursion: v-free
                self.ext_tgt[c] = tc
                self.ext_mul[c] = k % self.p
        return c

    def row_of_poly(self, poly, prefix=()):
        r"""
        Coefficient row of an encoded-ring polynomial, with an optional
        exponent ``prefix`` prepended to every exponent (used to encode
        ``U.ring`` polynomials as topforms via prepending ``u^n*v^0``;
        this would be accomplished by setting ``prefix = (n, 0)``).
        """
        col_get = self.col_of.get
        col_new = self.col
        row = {}
        if prefix:
            for e, c in poly.iterator_exp_coeff(as_ETuples=False):
                ee = prefix + e
                cc = col_get(ee)
                if cc is None:
                    cc = col_new(ee)
                row[cc] = int(c)
        else:
            for e, c in poly.iterator_exp_coeff(as_ETuples=False):
                cc = col_get(e)
                if cc is None:
                    cc = col_new(e)
                row[cc] = int(c)
        return row

    def row_deg(self, row):
        r"""Total degree of a row (`-1` for the zero row, as for polynomials)."""
        return max((self.deg[c] for c in row), default=-1)

    # ------------------------------------------------------------------
    # normal forms

    def ensure_nf(self, cols):
        r"""
        Make sure every column in ``cols`` has a known normal form.  New
        reducible monomials are reduced by Singular in one batch per round;
        a fix-point substitution pass afterwards guarantees the stored rows
        are supported on normal monomials only (i.e. they are the canonical
        reduced normal forms), independently of how much tail reduction
        Singular performed.
        """
        nf = self.nf
        exps = self.exps
        J = self._jaclm_arr
        pending = [c for c in cols if nf[c] is None]
        newly = []
        while pending:
            unknown = [c for c in sorted(set(pending)) if nf[c] is None]
            pending = []
            if not unknown:
                break
            # Vectorized staircase-membership test for the whole round.
            E = np.array([exps[c] for c in unknown], dtype=np.int64)
            if J is None:
                normal_mask = np.ones(len(unknown), dtype=bool)
            else:
                normal_mask = ~(J[None, :, :] <= E[:, None, :]).all(axis=2).any(axis=1)
            todo = []
            for c, isn in zip(unknown, normal_mask):
                if isn:
                    nf[c] = True
                else:
                    todo.append(c)
            if not todo:
                break
            ring = self.U.xring
            monos = [ring.monomial(*exps[c]) for c in todo]
            red = list(ff.reduce(ring.ideal(monos), self.U.jac))
            for c, q in zip(todo, red):
                row = self.row_of_poly(q)
                nf[c] = row
                newly.append(c)
                pending.extend(row)
        # Fix-point substitution on the rows computed in this call: replace
        # any reducible monomial remaining in a stored normal form by its
        # own normal form.  The support of the normal form of a monomial
        # lies strictly below it in the term order, so the recursion is
        # well founded; when Singular fully tail reduces (the default),
        # this pass never fires, and rows stored by earlier calls are
        # already flat by this same invariant.
        p = self.p

        def _flat(c):
            val = nf[c]
            bad = [k for k in val if nf[k] is not True]
            if not bad:
                return val
            row = dict(val)
            while bad:
                k = bad.pop()
                sub = _flat(k)
                cv = row.pop(k, 0)
                if cv:
                    for k2, v2 in sub.items():
                        nv = (row.get(k2, 0) + cv * v2) % p
                        if nv:
                            row[k2] = nv
                        else:
                            row.pop(k2, None)
                bad = [kk for kk in row if nf[kk] is not True]
            nf[c] = row
            return row

        for c in newly:
            if nf[c] is not True:
                _flat(c)

    def nf_row(self, row):
        r"""Canonical normal form of a row (callers batch ``ensure_nf`` first)."""
        p = self.p
        nf = self.nf
        out = {}
        get = out.get
        # Accumulate without intermediate reductions (values stay far below
        # 2^63 for the prime sizes used) and reduce once at the end.
        for c, v in row.items():
            f = nf[c]
            if f is True:
                out[c] = get(c, 0) + v
            else:
                for k2, v2 in f.items():
                    out[k2] = get(k2, 0) + v * v2
        return {k: r for k, v in out.items() if (r := v % p)}

    # ------------------------------------------------------------------
    # linear steps

    def ext_row(self, row):
        r"""Apply the encoded exterior derivative (a scaled column scatter)."""
        p = self.p
        tgt = self.ext_tgt
        mul = self.ext_mul
        out = {}
        get = out.get
        for c, v in row.items():
            t = tgt[c]
            if t >= 0:
                out[t] = get(t, 0) + v * mul[c]
        return {k: r for k, v in out.items() if (r := v % p)}

    def step_image(self, c):
        r"""Return the cached sparse image ``d(NF(c))`` of column ``c``."""
        cached = self.step_arr.get(c)
        if cached is not None:
            return cached

        normal_form = self.nf[c]
        terms = ((c, 1),) if normal_form is True else normal_form.items()
        image = {}
        p = self.p
        ext_tgt = self.ext_tgt
        ext_mul = self.ext_mul
        for source, coefficient in terms:
            target = ext_tgt[source]
            if target < 0:
                continue
            value = (image.get(target, 0)
                     + coefficient * ext_mul[source]) % p
            if value:
                image[target] = value
            else:
                image.pop(target, None)

        cached = (tuple(image), tuple(image.values()))
        self.step_arr[c] = cached
        return cached

    def elementary_reduction_step(self, rows):
        r"""
        Apply one batched Griffiths-Dwork step to sparse coefficient rows.

        For large batches, normal-form substitution and exterior
        differentiation are composed into the cached column map
        ``c -> d(NF(c))``.  One sparse matrix product then computes the
        complete step without materializing an intermediate normal form or
        performing a second derivative-scatter conversion.  Small batches
        retain the lower-overhead direct dictionary implementation.
        """
        if not rows:
            return []
        allcols = set()
        nnz = 0
        for row in rows:
            allcols.update(row)
            nnz += len(row)
        self.ensure_nf(allcols)
        p = self.p

        # For tiny batches the sparse-matrix round trip costs more than the
        # direct dictionary substitution; both compute the same map.
        # TODO: replace this with adaptive threshold (how?)
        if nnz < 400:
            return [self.ext_row(self.nf_row(row)) for row in rows]

        # Assemble the batch and composed operator in the stable global
        # column coordinates.  Keeping those coordinates avoids the Python
        # remapping cost of compact per-call source and target indices.
        size = len(self.exps)
        xr, xc, xv = [], [], []
        for i, row in enumerate(rows):
            xr.extend([i] * len(row))
            xc.extend(row.keys())
            xv.extend(row.values())
        X = _sp.csr_matrix(
            (xv, (xr, xc)),
            shape=(len(rows), size),
            dtype=np.int64,
        )

        tr, tc, tv = [], [], []
        for column in allcols:
            image_columns, image_values = self.step_image(column)
            tr.extend([column] * len(image_columns))
            tc.extend(image_columns)
            tv.extend(image_values)
        if not tv:
            return [{} for _ in rows]
        T = _sp.csr_matrix(
            (tv, (tr, tc)),
            shape=(size, size),
            dtype=np.int64,
        )

        # X @ T accumulates at most one product per active source in an
        # input row.  With residues below 2^26, at most 2^11 products fit in
        # int64.  Wider rows retain the existing 13-bit coefficient split.
        maxrow = max(len(row) for row in rows)
        if maxrow <= 2048:
            Z = X @ T
            Z.data %= p
        else:
            half = 1 << 13
            Tlo = T.copy()
            Tlo.data %= half
            Tlo.eliminate_zeros()
            Thi = T.copy()
            Thi.data >>= 13
            Thi.eliminate_zeros()
            Z = X @ Tlo
            Zh = X @ Thi
            Z.data %= p
            Zh.data = (Zh.data % p) * (half % p) % p
            Z = Z + Zh
            Z.data %= p
        Z.eliminate_zeros()

        out = []
        indptr, indices, data = Z.indptr, Z.indices, Z.data
        for i in range(len(rows)):
            start, stop = indptr[i], indptr[i + 1]
            out.append(dict(zip(
                indices[start:stop].tolist(),
                data[start:stop].tolist(),
            )))
        return out

    def echelonize(self, rows, want_pivots=False):
        r"""
        Reduced row echelon basis of the span of ``rows``.

        The active columns are sorted by the cached degrevlex key in
        decreasing order -- the same ordering ``echelonized_basis_poly``
        uses -- so pivots (and the output) coincide exactly with the
        polynomial implementation.  Returns a list of ``(row, pivot_column)``
        pairs, plus the pivot rows of the input when ``want_pivots`` is set.
        """
        if not rows:
            return ([], []) if want_pivots else []
        active = set()
        for row in rows:
            active.update(row)
        cols = sorted(active, key=lambda c: self.skey[c], reverse=True)
        # Sparse entries-dict construction (the input rows are sparse), with
        # the column renumbering done as one vectorized lookup.
        colmap = np.empty(len(self.exps), dtype=np.int64)
        colmap[np.fromiter(cols, dtype=np.int64, count=len(cols))] = \
            np.arange(len(cols), dtype=np.int64)
        xr, xc, xv = [], [], []
        for i, row in enumerate(rows):
            xr.extend([i] * len(row))
            xc.extend(row.keys())
            xv.extend(row.values())
        xj = colmap[np.array(xc, dtype=np.int64)].tolist()
        entries = dict(zip(zip(xr, xj), xv))
        mat = Matrix(self.U.ring.base_ring(), len(rows), len(cols), entries,
                     sparse=False)
        piv = mat.pivot_rows() if want_pivots else None
        mat.echelonize()
        out = []
        for mrow in mat.rows():
            nz = mrow.nonzero_positions()
            if not nz:
                break
            row = {cols[j]: int(mrow[j]) for j in nz}
            out.append((row, cols[nz[0]]))
        return (out, piv) if want_pivots else out

    def linear_normal_form(self, rows, basis):
        r"""
        Eliminate, in place, the pivot monomial of each element of the
        echelonized ``basis`` (a list of ``(row, pivot_column)`` pairs)
        from every row of ``rows`` -- the row version of
        ``linear_normal_form_p``.
        """
        p = self.p
        pivot_rows = {pc: brow for brow, pc in basis}
        for row in rows:
            # Because ``basis`` is in reduced row echelon form, a basis row
            # is zero in every other pivot column.  Hence the required
            # coefficients can be captured from the target before any
            # subtraction, and pivots absent from this sparse row need never
            # be inspected.
            reductions = []
            for pc, c in row.items():
                brow = pivot_rows.get(pc)
                if brow is not None and c:
                    reductions.append((brow, c))

            for brow, c in reductions:
                for k, v in brow.items():
                    nv = (row.get(k, 0) - c * v) % p
                    if nv:
                        row[k] = nv
                    else:
                        row.pop(k, None)

    # ------------------------------------------------------------------
    # syzygy slices and the U^r_q recursion

    def syz_rows(self, deg):
        r"""
        Coefficient rows of the degree ``deg`` slice of the non-trivial
        syzygy space, together with the leading monomial exponent of each
        lift.  Uses the same cached staircase enumeration and armed-profile
        filter as ``_basis_syzygies_with_lms``; a staircase lift is the
        exponent shift of a syzygy generator, so no polynomial products
        are formed.
        """
        from .rham_koszul import _syz_quotients
        U = self.U
        quotients, syzlm_exps = _syz_quotients(U, deg)
        # filt = None
        # if U.profile is not None and U.profile.get("armed") and U.r == 2: # TODO
        #     filt = U.profile["syzlm"]
        # numpy views of the syzygy term exponents, for vectorized shifts.
        syz_np = getattr(self, "_syz_np", None)
        if syz_np is None:
            syz_np = [(np.array([e for e, _ in terms], dtype=np.int64),
                       [c for _, c in terms])
                      for terms in self.syz_terms]
            self._syz_np = syz_np

        col_get = self.col_of.get
        col_new = self.col
        rows = []
        liftlms = []
        for qe, i in quotients:
            liftlm = tuple(q + s for q, s in zip(qe, syzlm_exps[i]))
            # if filt is not None and liftlm not in filt:
            #     continue
            T, C = syz_np[i]
            shifted = (T + np.array(qe, dtype=np.int64)).tolist()
            row = {}
            for ee, c in zip(map(tuple, shifted), C):
                cc = col_get(ee)
                if cc is None:
                    cc = col_new(ee)
                row[cc] = c
            rows.append(row)
            liftlms.append(liftlm)
        return rows, liftlms

    def compute_basis_U(self, r, q):
        r"""
        Row version of ``basis_U``: fill the engine caches for the degree
        ``q`` slice at reduction order ``r``, recording pivot provenance for
        the profile exactly as the polynomial implementation does.
        """
        if r < 0:
            raise ValueError("r must be non-negative")
        if (r, q) in self.basis_U:
            return
        if r == 0:
            self.basis_U[(r, q)] = []
            return
        U = self.U
        if r == 1:
            syz, liftlms = self.syz_rows(q - U.deg)
            self.basisWdown[(r, q)] = [self.ext_row(row) for row in syz]
            self.wdown_liftlms[(r, q)] = liftlms
            self.basis_U[(r, q)] = []
            return
        self.compute_basis_U(r - 1, q)
        self.compute_basis_U(r - 1, q + U.deg)
        rels = self.elementary_reduction_step(self.basisWdown[(r - 1, q + U.deg)])
        inputs = [row for row, _ in self.basis_U[(r - 1, q)]] + rels
        # if r == 2 and U.record_profile:
        #     new, piv = self.echelonize(inputs, want_pivots=True)
        #     off = len(self.basis_U[(r - 1, q)])
        #     liftlms = self.wdown_liftlms.get((r - 1, q + U.deg), [])
        #     for i in piv:
        #         if i >= off and i - off < len(liftlms):
        #             U.used_syzlms.add(liftlms[i - off])
        # else:
        new = self.echelonize(inputs) # TODO: When you uncomment the above, indent this so it's in the else.
        basis_q = []
        wdown = []
        for row, pc in new:
            # The pivot is the leading monomial in degree-compatible
            # degrevlex order, hence its degree is the degree of the row.
            pivot_degree = self.deg[pc]
            if pivot_degree == q:
                basis_q.append((row, pc))
            elif pivot_degree < q:
                wdown.append(row)
        self.basis_U[(r, q)] = basis_q
        self.basisWdown[(r, q)] = wdown

    def hom_reduce(self, rows, r):
        r"""
        Row version of ``_hom_reduce_helper``: repeatedly apply one batched
        Griffiths-Dwork step and, for ``r > 1``, eliminate against the
        precomputed slice $U^r_q$, lowering ``q`` by ``deg f`` each round.
        The list ``rows`` is rewritten in place.
        """
        if not rows:
            return
        q = max(self.row_deg(row) for row in rows)
        while q > 0:
            rows[:] = self.elementary_reduction_step(rows)
            # Zero is fixed by all remaining linear reduction steps.  In
            # particular, no lower-degree U-slices need to be constructed.
            if not any(rows):
                return
            if r > 1:
                self.compute_basis_U(r, q)
                self.linear_normal_form(rows, self.basis_U[(r, q)])
            q -= self.U.deg

    # ------------------------------------------------------------------
    # encode / decode

    def encode(self, poly):
        r"""Row of ``tox(poly) * u^n`` for a polynomial of ``U.ring``."""
        return self.row_of_poly(poly, prefix=(self.U.dim, 0))

    def decode(self, row):
        r"""
        Row version of ``fromx``: keep the ``v``-free part, drop the ``u``
        exponent, and merge coefficients on equal ``x``-parts.  Keys of the
        result are exponent tuples of ``U.ring``.
        """
        p = self.p
        out = {}
        for c, v in row.items():
            e = self.exps[c]
            if e[1]:
                continue
            k = e[2:]
            nv = (out.get(k, 0) + v) % p
            if nv:
                out[k] = nv
            else:
                out.pop(k, None)
        return out

    def xkey(self, e):
        r"""Degrevlex sort key of an ``U.ring`` exponent tuple (cached)."""
        k = self._xskey_cache.get(e)
        if k is None:
            k = self._xskey_fn(ETuple(list(e)))
            self._xskey_cache[e] = k
        return k

########## End of the LinearEngine class ##########


def engine_gauss_manin_helper(U, der, L, ret=None):
    r"""
    Drop-in replacement for the body of ``gauss_manin_helper`` over prime
    finite fields, using a ``LinearEngine`` attached to ``U``.  It performs
    the identical monomial closure -- same frontier ordering, same caps,
    same degree checks -- and returns the same ``RKGaussManinData`` record
    (basis monomials, exponent key, projection and action matrices).

    This is an internal function for sage_periods, and is not meant to be
    called by the user.
    """
    from .rham_koszul import RKGaussManinData

    eng = U.eng
    if eng is None:
        eng = LinearEngine(U)
        U.eng = eng
    n = U.dim
    F = U.ring.base_ring()
    p = eng.p
    r = U.r

    # Reduce the initial generators; keep the decoded rows for the
    # projection matrix.
    enc = [eng.encode(poly) for poly in L]
    eng.hom_reduce(enc, r)
    projx = [eng.decode(row) for row in enc]

    # Terms of -der, pre-encoded: multiplying the encoded monomial of a
    # basis element x^E by -der only shifts exponents.
    negder = [(tuple(e), (-int(c)) % p) for e, c in der.dict().items()]

    basis = []      # x-exponent tuples, in the exact legacy order
    basis_set = set()
    gmx = []        # decoded reduced images, one per basis element
    max_gm_degree = -1
    frontier = set()
    for row in projx:
        frontier.update(row)
    xmons = list(frontier)

    while xmons:
        xmons = sorted(xmons, key=eng.xkey)
        basis.extend(xmons)
        basis_set.update(xmons)

        # Images of the frontier under multiplication by -der, encoded
        # directly by exponent shifts, then reduced in batch.
        nf = []
        for Em in xmons:
            row = {}
            for Ed, c in negder:
                ee = (n, 0) + tuple(a + b for a, b in zip(Ed, Em))
                cc = eng.col(ee)
                row[cc] = c
            nf.append(row)
        eng.hom_reduce(nf, r)
        newx = [eng.decode(row) for row in nf]
        gmx.extend(newx)

        newmons = set()
        for row in newx:
            for exponent in row:
                newmons.add(exponent)
                degree = sum(exponent)
                if degree > max_gm_degree:
                    max_gm_degree = degree

        if max_gm_degree > (U.dim - 1) * U.deg:
            raise ReductionOrderTooSmallError("Linear Engine: The reduction order escaped the filtration bound.",r)

        xmons = list(newmons.difference(basis_set))

    ret = RKGaussManinData() if ret is None else ret
    ret.basis = [U.ring.monomial(*e) for e in basis]
    ret.ebasis = tuple((tuple(e),) for e in basis)

    nb = len(basis)
    xindex = {e: i for i, e in enumerate(basis)}
    gm_entries = {}
    for j, row in enumerate(gmx):
        for e, v in row.items():
            gm_entries[(xindex[e], j)] = v
    ret.gm = matrix(F, nb, nb, gm_entries)
    proj_entries = {}
    for j, row in enumerate(projx):
        for e, v in row.items():
            proj_entries[(xindex[e], j)] = v
    ret.proj = matrix(F, nb, len(L), proj_entries)
    return ret
