r"""Code for creation and verification of partial certificates, a la Section 7.3 of Lairez"""

from sage.all import *

from dataclasses import dataclass

@dataclass
class PeriodCertificate:
    r"""
    Container for a partial certificate in the sense of Lairez 2016,
    Section 7.3.

    Fields:

    * ``a``, ``f`` -- the prepared homogeneous numerator and squarefree
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

    def twisted_differential(encoded):
        """
        Take an encoded n-form and return the coefficient of D_f(form),
        a topform. Assumes a valid encoded n-form is passed in.
        """
        components = [A.zero()]*n

        for exponents, coefficient in encoded.items():
            e = tuple(exponents)  # (u exponent, v exponent, x exponents )
            components[e[1] - 1] += A({e[2:]: K(coefficient)})

        return sum((b.derivative(x) - b * cert.f.derivative(x)
                    for b, x in zip(components, A.gens())), A.zero())

    try:
        # build M, verify B, rho, gammas are right shape
        M = []
        for exponents in cert.basis_exponents:
            e = tuple(exponents)
            if len(e) != n or (sum(e) + n) % N:
                return False
            M.append(A({e: K.one()}))
        Msize = len(M)
        if (len(cert.rho0) != Msize or cert.B.dimensions() != (Msize, Msize)
                or len(cert.gammas) != Msize):
            return False

        rho0 = vector(K, cert.rho0)
        rho0_poly = sum((rho0[j] * M[j] for j in range(Msize)), A.zero())

        # CHECK 1: Initial reduction: rho0 = a + D_f(beta0).
        if rho0_poly != cert.a + twisted_differential(cert.beta0):
            return False

        # CHECK 2: Each column reduces -f_delta times its basis elt of M.
        f_delta = A({e: c.derivative() for e, c in cert.f.dict().items()})
        for j, mu in enumerate(M):
            column = sum((K(cert.B[i, j]) * M[i]
                          for i in range(Msize)), A.zero())
            if column != -f_delta * mu + twisted_differential(cert.gammas[j]):
                return False

        coeffs = [K(c) for c in cert.operator_coefficients]
        if not coeffs or not any(coeffs):
            return False

        # CHECK 3: linear dependency among the operator coefficients
        if Msize:
            relation = vector(K, Msize)
            rho = rho0
            for c in coeffs:
                relation += c * rho
                rho = vector(K, [entry.derivative() for entry in rho]) + cert.B * rho
            if not relation.is_zero():
                return False
        return True
    except (ArithmeticError, TypeError, ValueError, IndexError, AttributeError):
        return False


__all__ = ["PeriodCertificate", "verify_certificate"]
