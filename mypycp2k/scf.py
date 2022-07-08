def add_ot_never_fail(scf):
    """
    from Jorge. Very stable, very slow
    :param scf:
    :return:
    """
    OT = scf.OT
    OT.Minimizer = 'CG'
    OT.Preconditioner = 'FULL_KINETIC'
    OT.Algorithm = 'STRICT'
    OT.Ortho_irac = 'POLY'

def add_diagonalization(scf):
    DIAGONALIZATION = scf.DIAGONALIZATION
    DIAGONALIZATION.Section_parameters = 'T'


def add_smear(scf, method='FERMI_DIRAC', electronic_temperature=300):
    SMEAR = scf.SMEAR
    SMEAR.Method = method
    SMEAR.Electronic_temperature = electronic_temperature


def add_mixing(scf, method='BROYDEN_MIXING'):
    MIXING = scf.MIXING
    MIXING.Method = method
    MIXING.Alpha = 0.2
    MIXING.Beta = 1.5
    MIXING.Nbroyden = 8


def add_mos(scf, added_mos=100):
    scf.Added_mos = added_mos


def add_broyden_mixing(scf):
    """
    diagonalization with mixing and smear
    settings from here: https://www.cp2k.org/exercises:2017_uzh_cmest:pdos
    :param scf:
    :return:
    """
    add_diagonalization(scf)
    scf.DIAGONALIZATION.Eps_adapt = 0.01
    scf.DIAGONALIZATION.Algorithm = 'STANDARD'
    add_smear(scf)
    add_mixing(scf)
    add_mos(scf)


def add_ot(scf,
           minimizer='DIIS',
           preconditioner='FULL_ALL',
           energy_gap=0.001,
           stepsize=0.05):
    """
    seems to perform best for molecules so far
    :param scf:
    :param minimizer:
    :param preconditioner:
    :param energy_gap:
    :param stepsize:
    :return:
    """
    OT = scf.OT
    OT.Minimizer = minimizer
    OT.Preconditioner = preconditioner  # FULL_SINGLE_INVERSE
    OT.Energy_gap = energy_gap
    OT.Stepsize = stepsize  # best for mol. 100100


def remove_ot(scf):
    """
    set all OT nonprotected attr to zero
    """

    OT = scf.OT
    OT_dict = OT.__dict__
    all_keys = list(OT_dict.keys())
    non_protected_keys = [key for key in all_keys if not key.startswith('_')]

    for key in non_protected_keys:
        setattr(OT, key, None)


def put_all_dict_values_to_zero():
    pass
