r"""Code for creation and verification of partial certificates, a la Section 7.3 of Lairez"""

from sage.all import *

from dataclasses import dataclass

@dataclass
class PeriodCertificate:
    r"""
    Container for a partial certificate in the sense of Lairez 2016,
    Section 7.3.

    Fields:

    * ``a``, ``f`` -- the prepared homogeneous numerator and
      denominator, as polynomials over $\QQ(t)$;
    * ``basis_exponents`` -- the exponent vectors of the monomial basis
      $\mathcal M$;
    * ``rho0`` -- the coordinates of the reduction of $a\omega$ in
      $\mathcal M$, over $\QQ(t)$;
    * ``B`` -- the matrix over $\QQ(t)$ of the map
      $\rho \mapsto [-f^{\delta}\rho]$ in the basis $\mathcal M$;
    * ``beta0`` -- the encoded cofactor of the reduction of $a\omega$, as a
      dictionary mapping exponent vectors of the auxiliary $u,v$-ring to
      elements of $\QQ(t)$;
    * ``gammas`` -- one such dictionary per basis monomial, the cofactor of
      the reduction of $-f^{\delta}\mu\,\omega$;
    * ``operator_coefficients`` -- the dense coefficient list
      $[c_0, \dots, c_m]$ of the certified operator
      $L = \sum_k c_k \delta^k$;
    """
    a: object
    f: object
    basis_exponents: tuple
    rho0: list
    B: object
    beta0: dict
    gammas: list
    operator_coefficients: list | None


def verify_certificate(cert):
    r"""
    Verify a ``PeriodCertificate`` ``cert`` by exact computation over $\QQ(t)$.

    Check the initial reduction, each Gauss-Manin matrix column, and the
    operator relation by exact arithmetic over $\QQ(t)$.
    
    A successful verification proves that the operator formed from
    ``cert.operator_coefficients`` annihilates all the periods of the
    prepared fraction $a/f^{q}$. The proof does not rely on the algorithm
    that produced the certificate. This routine checks the coefficients in
    ``cert.operator_coefficients``; any separately minimized operator needs
    its own verification.

    INPUT:

    * ``cert`` -- A ``PeriodCertificate``.

    OUTPUT:

    * A boolean.

    EXAMPLES:

    A hand-built certificate for the constant period of $1/(x_0x_1)$:

        sage: K.<t> = FractionField(PolynomialRing(QQ, 't'))
        sage: A.<x0, x1> = PolynomialRing(K)
        sage: cert = PeriodCertificate(x0*x1, x0*x1, ((0, 0),), [1],
        ....:     matrix(K, [[0]]), {(1, 1, 1, 0): 1}, [{}], [0, 1])
        sage: verify_certificate(cert)
        True
    """
    A = cert.f.parent()
    K = A.base_ring()
    n = A.ngens()
    N = cert.f.degree()
    # cheap, dumb tests
    if (not cert.f or N <= 0 or not cert.f.is_homogeneous()
            or cert.a.parent() != A
            or (cert.a and (not cert.a.is_homogeneous()
                            or (cert.a.degree() + n) % N))):
        return False

    # Instead of running checks over A = QQ(t)[x_0,...,x_n], multiply all polys
    # by an lcm and run over S = QQ[x_0, ..., x_n, t].
    # This avoids many Singular gcd calls.
    R = K.ring()
    S = PolynomialRing(K.base_ring(), list(A.variable_names()) + [R.variable_name()])
    to_S = R.hom([S.gen(n)], S)

    def denominator_lcm(coefficients):
        return lcm([K(c).denominator() for c in coefficients if c] or [R.one()])

    def cleared(terms, denom):
        """``denom * sum(c * x^e)`` in S; ``denom*c`` must be a polynomial in t."""
        out = {}
        for e, c in terms:
            q = denom * K(c)
            if q.denominator() != 1:
                raise ArithmeticError("denominator does not divide the multiplier")
            for k, ck in enumerate(q.numerator().list()):
                if ck:
                    out[tuple(e) + (k,)] = ck
        return S(out)

    # cast f, its denominator, and f^\delta into S.
    fdict = cert.f.dict()
    denom_f = denominator_lcm(fdict.values())
    f_S = cleared(fdict.items(), denom_f)
    df_S = [f_S.derivative(S.gen(i)) for i in range(n)]
    denom_f_S = to_S(denom_f)

    def check_reduction(lhs, encoded):
        """Exact test of ``lhs == D_f(encoded)``; ``lhs`` is a dict {exponent : QQ(t)}."""
        components = [[] for _ in range(n)]
        for exponents, coefficient in encoded.items():
            if coefficient:
                e = tuple(exponents)  # (u exponent, v exponent, x exponents...)
                components[e[1] - 1].append((e[2:], coefficient))
        denom_b = denominator_lcm(encoded.values())
        # denom_b*denom_f*D_f(beta) = sum_i (denom_f d_i b~_i - b~_i d_i f~), with b~ = denom_b*b, f~ = denom_f*f.
        # Build rhs
        rhs = S.zero()
        for i, terms in enumerate(components):
            if terms:
                b = cleared(terms, denom_b)
                rhs += denom_f_S * b.derivative(S.gen(i)) - b * df_S[i]
        # Build lhs denominator
        D = lcm(denom_b * denom_f, denominator_lcm(lhs.values()))
        return cleared(lhs.items(), D) == to_S(D // (denom_b * denom_f)) * rhs

    try:
        # Build M
        M = []
        for exponents in cert.basis_exponents:
            e = tuple(exponents)
            if len(e) != n or (sum(e) + n) % N:
                return False
            M.append(e)
        Msize = len(M)
        # Verify B, rho, gammas are right shape
        if (len(cert.rho0) != Msize or cert.B.dimensions() != (Msize, Msize)
                or len(cert.gammas) != Msize):
            return False

        # CHECK 1: Initial reduction: rho0 = a + D_f(beta0).
        lhs = {}
        for e, c in zip(M, cert.rho0):
            lhs[e] = lhs.get(e, 0) + K(c)
        for e, c in cert.a.dict().items():
            e = tuple(e)  # ETuple keys do not match tuple keys
            lhs[e] = lhs.get(e, 0) - c
        if not check_reduction(lhs, cert.beta0):
            return False

        # CHECK 2: Each column reduces -f_delta times its basis elt of M.
        f_delta = [(tuple(e), c.derivative()) for e, c in fdict.items()]
        for j, mu in enumerate(M):
            lhs = {}
            for i, e in enumerate(M):
                lhs[e] = lhs.get(e, 0) + K(cert.B[i, j])
            for e, c in f_delta:
                emu = tuple(x + y for x, y in zip(e, mu))
                lhs[emu] = lhs.get(emu, 0) + c
            if not check_reduction(lhs, cert.gammas[j]):
                return False

        coeffs = [K(c) for c in cert.operator_coefficients]
        if not coeffs or not any(coeffs):
            return False

        # CHECK 3: linear dependency among the operator coefficients.
        # Take the lcm g of all denominators of rho0, B and the c_k first.
        if Msize:
            g = denominator_lcm(list(cert.rho0) + cert.B.list() + coeffs)
            G = matrix(R, Msize, Msize, [R(g * K(c)) for c in cert.B.list()])
            r = vector(R, [R(g * K(c)) for c in cert.rho0])
            c_R = [R(g * c) for c in coeffs]
            m = len(coeffs) - 1
            dg = g.derivative()
            relation = vector(R, Msize)
            for k, c in enumerate(c_R):
                relation += (c * g**(m - k)) * r
                if k < m:
                    r = vector(R, [g * x.derivative() - (k + 1) * dg * x for x in r]) + G * r
            if not relation.is_zero():
                return False
        return True
    except (ArithmeticError, TypeError, ValueError, IndexError, AttributeError):
        return False


__all__ = ["PeriodCertificate", "verify_certificate"]
