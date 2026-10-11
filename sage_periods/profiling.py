r"""(Generalized) Syzygy profiling, compatible with both polynomial and linear engines."""
from dataclasses import dataclass


@dataclass(slots=True)
class ProfileData:
    r"""
    Container storing the syzygy and trivial syzygy data shared across primes
    and evaluation points, used to speed up computation of the Gauss-Manin
    connection at reduction order ``r = 2``.

    A single instance is created in ``compute_period_annihilator`` and is read
    and modified in place by the probe, ``compute_gauss_manin_connection`` and
    ``RhamKoszulData``.  All stored monomial data are exponent tuples, so the
    profile is independent of the coefficient field and can be shared across
    primes.

    States:

    * ``"disabled"`` -- profiling is not in use (initial state, and ``r != 2``).
    * ``"recording"`` -- each successful point produces a candidate profile;
      no filtering is applied.
    * ``"armed"`` -- two consecutive candidates agreed; the accepted profile is
      frozen and used to skip the trivial-syzygy Groebner basis and to filter
      the syzygy staircase lifts.

    A candidate profile is the tuple ``(tsyzlms, syzgen_lms, syzlms)``:

    * ``tsyzlms`` -- leading monomial exponents of the trivial-syzygy Groebner basis;
    * ``syzgen_lms`` -- leading monomial exponents of the syzygy generators ``U.syz``;
    * ``syzlms`` -- ``frozenset`` of leading monomial exponents of the staircase
      lifts which contributed a pivot.

    This class should be compatible with both the polynomial-based reduction
    algorithm used in ``rham_koszul.py`` and the linear coefficient-based
    version of ``linear_closure.py``.
    """
    state: str = "disabled"
    # Most recent candidate; we arm when the next candidate agrees with it.
    prev_prof: tuple | None = None
    # Prime at which the accepted profile was produced. If this prime is later
    # discarded (e.g. in Stage 3), the profile must be reset.
    prof_prime: object = None
    # Accepted profile.
    tsyzlms: tuple = ()
    syzgen_lms: tuple = ()
    syzlms: frozenset = frozenset()

    # TODO: If this version is slow, instead of checking the entire profile on every candidate add,
    # ...just test a certain randomly generated small subset of the entries?

    def reset(self):
        r"""Reset ``ProfileData`` object to its initial state, in place."""
        self.state = "disabled"
        self.prev_prof = None
        self.prof_prime = None
        self.tsyzlms = ()
        self.syzgen_lms = ()
        self.syzlms = frozenset()

    def copy_from(self, other):
        r"""Overwrite this profile, in place, with the contents of ``other``."""
        self.state = other.state
        self.prev_prof = other.prev_prof
        self.prof_prime = other.prof_prime
        self.tsyzlms = other.tsyzlms
        self.syzgen_lms = other.syzgen_lms
        self.syzlms = other.syzlms

    def add_profile_candidate(self, profile, prime):
        r"""
        Stability / profile-selection manager for candidate profiles. Takes a new
        ``profile`` candidate ``(tsyzlms, syzgen_lms, syzlms)`` computed at one
        successful evaluation point modulo ``prime``, and arms the profile when
        it agrees with the previous candidate.

        Similar in spirit to ``recon_add_rat``.
        """
        if self.state != "recording":
            return

        # NOTE: Would insert first-candidate random sample generation here,
        # ...if we have to go that route.

        if self.prev_prof is not None and self.prev_prof == profile:
            self.tsyzlms, self.syzgen_lms, self.syzlms = profile
            self.prof_prime = prime
            self.state = "armed"
        else:
            self.prev_prof = profile
