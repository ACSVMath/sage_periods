r"""(Generalized) Syzygy profiling, compatible with both polynomial and linear engines."""
from dataclasses import dataclass, field

@dataclass(slots=True)
class ProfileData:
    r"""
    Container storing the syzygy and trivial syzygy data across primes & evaluation points,
    used to speed up computation of Gauss-Manin connection.

    This record keeps track of candidate profiles, the state-tracking information
    for our profiler, as well as the accepted candidate profile itself.
    This class should be compatible with both the polynomial-based reduction 
    algorithm used in ``rham_koszul.py`` and the linear coefficient-based 
    version of ``linear_closure.py``
    
    , compressed test values, current rational reconstructions,
    the set of already used points, and a candidate once reconstruction has stabilized.
    """
    state: str = field(default_factory=str) # Either "disabled", "recording", or "armed"
    # candidates: dict = field(default_factory=dict) # Keyed by either points or primes; 
    # values are either dicts of polynomial monomials, in the polynomial pipeline case, 
    # or dicts of monomial columns, in the linear pipeline case.
    prev_prof: tuple=field(default_factory=tuple) # Even simpler: instead of sorting candidates into "buckets", we move forward if we simply
    # agree with the last profile passed in.
    tsyzlms: dict=field(default_factory=dict) # Our final (accepted) candidate, trivial syzygies
    syzlms: dict=field(default_factory=dict) # Our final (accepted) candidate, all syzygies

    # TODO: If this version is slow, instead of checking the entire profile on every candidate add,
    # ...just test a certain randomly generated small subset of the entries?

    def __init__(self):         # TODO: build this out?
        r"""Initialize ``ProfileData`` to its default state."""
        self.state = "disabled"
        # self.candidates = {}
        self.prev_prof = {}
        self.tsyzlms = {}
        self.syzlms = {}
        self.first_add = False
    
    def reset(self):
        r"""Rest ``ProfileData`` object to its initial state"""
        self.state = "disabled"
        self.candidates = {}
        self.tsyzlms = {}
        self.syzlms = {}
    
    def add_profile_candidate(self,profile):
        r"""
        Stability / profile-selection manager for candidate profiles. Takes a new
        ``profile`` object (tsyzlms, syzlms) -- (literally a tuple!) -- and 
        runs stability testing to see if it has occurred before,
        keeping that profile when it has.
        
        Similar in spirit to ``recon_add_rat``.
        """

        # NOTE: Would insert first-candidate random sample generation here,
        # ...if we have to go that route.

        if not self.prev_prof:
            self.prev_prof = profile
        else:
            if self.prev_profile == profile:
                self.tsyzlms = profile[0]
                self.syzlms = profile[1]
            
            self.prev_prof = profile

        return

