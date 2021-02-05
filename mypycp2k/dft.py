def set_dft(DFT,
            potential_file_name,
            basis_set_file_name,
            multiplicity=1,
            charge=0,
            uks=False):
    # del DFT.Basis_set_file_name
    DFT.Potential_file_name = potential_file_name
    DFT.Basis_set_file_name = basis_set_file_name
    DFT.Multiplicity = multiplicity
    DFT.Uks = uks
    DFT.Charge = charge


def set_scf(DFT, eps_scf=1.0E-6, max_scf=100):
    SCF = DFT.SCF
    SCF.Scf_guess = 'RESTART'  # atomic will be used if not exist
    SCF.Eps_scf = eps_scf
    SCF.Max_scf = max_scf

def set_nonperiodic_poisson(DFT, solver='MT'):
    POISSON = DFT.POISSON
    POISSON.Periodic = 'NONE'
    POISSON.Poisson_solver = solver


def set_cutoff(DFT, cutoff=900, rel_cutoff=60, ngrids=5):
    MGGRID = DFT.MGRID
    MGGRID.Cutoff = cutoff
    MGGRID.Rel_cutoff = rel_cutoff
    MGGRID.Ngrids = ngrids

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


def set_qs(DFT,
           method='GAPW',
           eps_default=1.0E-15,
           eps_pgf_orb=1.E-200):
    QS = DFT.QS
    QS.Method = method
    QS.Eps_default = eps_default  # default: 1E-10. Need higher???
    QS.Eps_pgf_orb = eps_pgf_orb  # Overrides SQRT(EPS_DEFAULT) value. No default value. todo: clarify
