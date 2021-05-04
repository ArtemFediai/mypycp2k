class Error(Exception):
    """Base class for other exceptions"""
    pass


class SCQPSolutionNotFound(Error):
    """Raised when Self-consistent quasi-particle solution not found (GW)"""

    def __init__(self, message='Self-consistent quasi-particle solution not found (GW)'):
        # self.message = message
        super().__init__(message)


class SCFNotConvergedNotPossibleToRunMP2(Error):
    """SCF not converged: not possible to run MP2"""

    def __init__(self, message='SCF not converged: not possible to run MP2'):
        # self.message = message
        super().__init__(message)
