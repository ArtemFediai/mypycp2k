def set_dft(DFT,
            potential_file_name,
            basis_set_file_name,
            multiplicity=1,
            charge=0,
            uks=False):
    DFT.Potential_file_name = potential_file_name
    DFT.Basis_set_file_name = basis_set_file_name
    DFT.Multiplicity = multiplicity
    DFT.Uks = uks
    DFT.Charge = charge


def set_basis_set_file_name(DFT,
                            basis_set_file_name):
    # del DFT.Basis_set_file_name
    DFT.Basis_set_file_name = basis_set_file_name


def set_scf(DFT, eps_scf=1.0E-6, max_scf=100, scf_guess='RESTART'):
    SCF = DFT.SCF
    SCF.Scf_guess = scf_guess  # atomic will be used if not exist
    SCF.Eps_scf = eps_scf
    SCF.Max_scf = max_scf


def set_nonperiodic_poisson(DFT, solver='MT'):
    POISSON = DFT.POISSON
    POISSON.Periodic = 'NONE'
    POISSON.Poisson_solver = solver


def set_periodic_poisson(DFT, solver='PERIODIC'):
    POISSON = DFT.POISSON
    POISSON.Periodic = 'XYZ'
    POISSON.Poisson_solver = solver


def set_cutoff(DFT, cutoff=900, rel_cutoff=60, ngrids=5, commensurate='.FALSE.'):
    MGRID = DFT.MGRID
    MGRID.Cutoff = cutoff
    MGRID.Rel_cutoff = rel_cutoff
    MGRID.Ngrids = ngrids
    MGRID.Commensurate = commensurate


def print_mo_cubes(PRINT, nhomo=5, nlumo=5, filename='=HOMO.txt'):
    MO_CUBES = PRINT.MO_CUBES
    MO_CUBES.Nhomo = nhomo
    MO_CUBES.Nlumo = nlumo
    MO_CUBES.Write_cube = False
    MO_CUBES.Filename = filename


def print_mo(PRINT, filename='MO', range='1 10'):
    MO = PRINT.MO
    MO.Eigenvectors = True
    MO.Filename = filename
    MO.Range = range


# v9 means for qp. or for cp2k ver.9.1.
def set_qs_v9(
        DFT,
        method='GAPW',
        eps_default=1.0E-10,
        eps_pgf_orb=None,
        eps_ppl=1E-2,
        epssvd=1E-8,
        epsrho0=1E-6,
        epsiso=1E-12,
        epsfit=1E-4
):
    QS = DFT.QS
    QS.Method = method
    QS.Eps_default = eps_default  # default: 1E-10. Need higher???
    QS.Eps_ppl = eps_ppl
    QS.Epsrho0 = epsrho0
    QS.Epssvd = epssvd
    QS.Epsiso = epsiso
    QS.Epsfit = epsfit
    if eps_pgf_orb is None:
        pass  # this will set it to sqrt(eps_default)
    else:
        QS.Eps_pgf_orb = eps_pgf_orb  # Overrides SQRT(EPS_DEFAULT) value. No default value. todo: clarify
    if method == 'LRIGPW':
        LRIGPW = QS.LRIGPW
        LRIGPW.Lri_overlap_matrix = 'INVERSE'
        LRIGPW.Shg_lri_integrals = '.TRUE.'
        # LRIGPW.Ppl_ri = '.TRUE.'  # Use LRI/RI for local pseudopotential
        # LRIGPW.Eps_o3_int = 1.00000000E-20  # threashold. may be sensitive


def set_qs(DFT,
           method='GAPW',
           eps_default=1.0E-10,
           eps_pgf_orb=1.E-200):
    QS = DFT.QS
    QS.Method = method
    QS.Eps_default = eps_default  # default: 1E-10. Need higher???
    QS.Eps_pgf_orb = eps_pgf_orb  # Overrides SQRT(EPS_DEFAULT) value. No default value. todo: clarify
