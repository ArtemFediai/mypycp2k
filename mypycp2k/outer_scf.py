def add_outer_scf(outer_scf, eps_scf=1.0E-10, max_scf=20):
    outer_scf.Eps_scf = eps_scf
    outer_scf.Max_scf = max_scf
