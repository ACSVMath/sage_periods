r"""Code for creation and verification of partial certificates, a la Section 7.3 of Lairez"""

from sage.all import *

from dataclasses import dataclass, field

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
    * ``returned_operator`` -- the operator object returned to the user
      (which may be a certified-by-division minimization of $L$).
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

    Returns ``True`` if and only if an equivalent (encoded $\QQ(t)$-polynomial 
    differential form) version of the criteria in Section 7.3 of Lairez are met.

    TODO: Also check homogeneity?
    
    A successful verification proves that the operator formed from
    ``cert.operator_coefficients`` annihilates all the periods of the
    prepared fraction $a/f^{q}$, hence of the original input; the proof does
    not rely on any property of the algorithm that produced the certificate.

    INPUT:

    * ``cert`` -- A ``PeriodCertificate``.

    OUTPUT:

    * A boolean.

    EXAMPLES:

    A certificate produced by the pipeline verifies:

            sage: from sage_periods import compute_diagonal_annihilator, verify_certificate
            sage: var('x y')
            sage: L, cert = compute_diagonal_annihilator(1/(1-x-y), certificate=True)
            sage: verify_certificate(cert)
            True

    Tampering with any part of the certificate makes it fail:

            sage: cert.B[0,0] += 1
            sage: verify_certificate(cert)
            False
    """

    # STEP 0: Sanity checks: 
    # - TODO: The parent rings of L and cert agree; L and cert have same "shape"?
    # - Homogeneity and degree bookkeeping of the prepared fraction?
    # - elements of cert.rho0, M, B have correct grading / degree?

    # STEP 1: Verify rho_0 = a*omega + D_f(beta_0).

    # STEP 2: Verify, for each mu basis elt of M, that the column of B 
    # ...corresponding to mu equals -f^delta * mu * omega + D_f(eta_mu).

    # STEP 3: Verify that sum_k c_k rho_k = 0 with rho_{k+1} = rho_k' + B rho_k.
    

    pass

# ^^NOTE: This relies on the fact that these three conditions are equivalent
# ...to the ones given in the paper. Perhaps a proof of this equivalence 
# ...should be made explicitly somewhere? But it's not that difficult anyway.


__all__ = ["PeriodCertificate", "verify_certificate"]