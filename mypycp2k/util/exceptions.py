class Error(Exception):
    """Base class for other exceptions"""
    pass


class SCQPSolutionNotFound(Error):
    """Raised when Self-consistent quasi-particle solution not found (GW)"""

    def __init__(self, message='Self-consistent quasi-particle solution not found (GW)'):
        super().__init__(message)


class SCFNotConvergedNotPossibleToRunMP2(Error):
    """SCF not converged: not possible to run MP2"""

    def __init__(self, message='SCF not converged: not possible to run MP2'):
        super().__init__(message)


class NaNInGW(Error):
    """NaN in GW table (last frame)"""

    def __init__(self, message='SCF GW is not extracted: NaN in GW energies'):
        super().__init__(message)

class LargeSigc(Error):
    """Sigc is too large (last frame)"""

    def __init__(self, message='SCF GW is extracted BUT SEEMS TO BE WRONG: Unphysicaly large Sigc'):
        super().__init__(message)

class IterationLimit(Error):
    """20 iteration limit is reached. Not actually converged"""

    def __init__(self, message='SCF GW is extracted BUT SEEMS TO BE WRONG: Iteration limit (20) is reached'):
        super().__init__(message)
