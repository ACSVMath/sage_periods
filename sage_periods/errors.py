r"""Handling of errors occurring /passed across multiple files (rham_koszul, pipeline, preparation, picard_fuchs)"""


# NOTE:
# Not every expected outcome becomes an exception. In particular, points being "bad" will be handled within reconstruction.py
# In general, 
# - If the current loop can handle the outcome locally, return a status/result.
# - If an entire nested attempt must be abandoned, 
#         - raise a typed exception and 
#         - catch it at the loop owning that retry policy.


class SagePeriodsError(Exception):
    r"""Base class for package computation errors.

    TESTS:

        sage: from sage_periods.errors import SagePeriodsError
        sage: str(SagePeriodsError('failed'))
        'failed'
    """
    pass

# Typical usage:
# raise BadPrimeError(message="Prime fails normalization checks in CRT", prime=7, stage="lifting") from exc
class BadPrimeError(SagePeriodsError):
    r"""Signal that the current prime cannot be used for the computation.

    The ``prime`` and ``stage`` attributes identify the failed attempt.

    TESTS:

        sage: from sage_periods.errors import BadPrimeError
        sage: e = BadPrimeError('bad prime', 7, 'lifting')
        sage: (e.prime, e.stage)
        (7, 'lifting')
    """
    def __init__(self, message, prime, stage): # stage: "gauss_manin", "dependency", "lifting"
        r"""Create error for ``prime`` at ``stage``.

        INPUT:

        * ``message`` -- Explanation of the failure.
        * ``prime`` -- Rejected prime.
        * ``stage`` -- Stage name, e.g. ``"lifting"``.

        TESTS:

            sage: from sage_periods.errors import BadPrimeError
            sage: BadPrimeError('bad', 7, 'lifting').prime
            7
        """
        super().__init__(message)
        self.prime = prime
        self.stage = stage
    pass

class NeedMorePrimesError(SagePeriodsError):
    r"""Signal that we either:
        - pruned out too many primes in CRT, or 
        - need more to rationally reconstruct. 
    So we should get more primes.

    TESTS:

        sage: from sage_periods.errors import NeedMorePrimesError, SagePeriodsError
        sage: isinstance(NeedMorePrimesError('more primes'), SagePeriodsError)
        True
    """
    def __init__(self, message):
        r"""Create error with explanatory ``message``.

        TESTS:

            sage: from sage_periods.errors import NeedMorePrimesError
            sage: str(NeedMorePrimesError('more primes'))
            'more primes'
        """
        super().__init__(message)
    pass

class ReductionOrderTooSmallError(SagePeriodsError):
    r"""Signal that reduction order ``r`` is too small.

    Attempted order is available as ``r``.

    TESTS:

        sage: from sage_periods.errors import ReductionOrderTooSmallError
        sage: ReductionOrderTooSmallError('increase order', 2).r
        2
    """
    def __init__(self, message, r):
        r"""Create error for attempted reduction order ``r``.

        INPUT:

        * ``message`` -- Explanation of the failure.
        * ``r`` -- Attempted order.

        TESTS:

            sage: from sage_periods.errors import ReductionOrderTooSmallError
            sage: str(ReductionOrderTooSmallError('increase order', 2))
            'increase order'
        """
        super().__init__(message)
        self.r = r
    pass

