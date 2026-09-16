r"""Handling of errors occurring /passed across multiple files (rham_koszul, pipeline, preparation, picard_fuchs)"""


# NOTE:
# Not every expected outcome becomes an exception. In particular, points being "bad" will be handled within reconstruction.py
# In general, 
# - If the current loop can handle the outcome locally, return a status/result.
# - If an entire nested attempt must be abandoned, 
#         - raise a typed exception and 
#         - catch it at the loop owning that retry policy.


class SagePeriodsError(Exception):
    """Base class for package computation errors."""
    pass

# Typical usage:
# raise BadPrimeError(message="Prime fails normalization checks in CRT", prime=7, stage="lifting") from exc
class BadPrimeError(SagePeriodsError):
    """The current prime cannot be used."""
    def __init__(self, message, prime, stage): # stage: "gauss_manin", "dependency", "lifting"
        super().__init__(message)
        self.prime = prime
        self.stage = stage
    pass

class NeedMorePrimesError(SagePeriodsError):
    """Either we pruned out too many primes in CRT, or need more to rationally reconstruct. Get more primes."""
    def __init__(self, message):
        super().__init__(message)
    pass

class ReductionOrderTooSmallError(SagePeriodsError):
    """The reduction escaped the expected filtration bound."""
    def __init__(self, message, r):
        super().__init__(message)
        self.r = r
    pass

class ProbeBasisCapExceededError(SagePeriodsError):
    """A probe exceeded its permitted basis size."""
    # TODO
    pass

# Other possiblities: failed certification?