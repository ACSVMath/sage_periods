r"""Result record shared by the polynomial and linear closure engines."""

from dataclasses import dataclass, field


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

