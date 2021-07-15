def set_pbe(XC):
    """
    regular PBE functional
    :param XC:
    :return:
    """
    #### XC_FUNCTIONAL
    XC_FUNCTIONAL = XC.XC_FUNCTIONAL
    XC.XC_FUNCTIONAL.Section_parameters = 'PBE'


def set_pbe0(XC, max_memory=1500):
    """
    add sections xc.functional and xc.HF in accordance to PBE0
    :param XC: pycp2k object
    :param max_memory: float
    :return: changes XC
    """
    #### XC_FUNCTIONAL
    XC_FUNCTIONAL = XC.XC_FUNCTIONAL
    XC.XC_FUNCTIONAL.Section_parameters = 'PBE'
    XC_FUNCTIONAL.PBE.Scale_x = 0.75
    XC_FUNCTIONAL.PBE.Scale_c = 1.00

    #### HF ####
    HF = XC.HF_add()
    HF.Fraction = 0.25
    HF.SCREENING.Eps_schwarz = 1.0E-11
    HF.SCREENING.Screen_on_initial_p = False
    HF.MEMORY.Max_memory = max_memory
    HF.MEMORY.Eps_storage_scaling = 0.1


def add_vdw(XC,
            vdw_parameters_file='dftd3.dat',
            my_vdw_potential='PAIR_POTENTIAL',
            pot_type='DFTD3(BJ)',
            ref_fun='PBE0'):
    VDW_POTENTIAL = XC.VDW_POTENTIAL
    VDW_POTENTIAL.Potential_type = my_vdw_potential
    PAIR_POTENTIAL1 = VDW_POTENTIAL.PAIR_POTENTIAL_add()
    PAIR_POTENTIAL1.Reference_functional = ref_fun  # pbd0?
    PAIR_POTENTIAL1.Parameter_file_name = vdw_parameters_file
    PAIR_POTENTIAL1.Type = pot_type
    PAIR_POTENTIAL1.PRINT_DFTD.Filename = 'vdw_out.dftd'


def add_gw_ver_0(xc,
                 wf_corr_num_proc=1,
                 rpa_num_quad_points=100,
                 size_freq_integ_group=-1,
                 corr_occ=10,
                 corr_virt=10,
                 ev_sc_iter=1,
                 max_memory_hf=500,
                 max_memory_wf=2000,
                 eps_schwarz = 1.0E-11):
    """
    my first version of GW settings.
    :param xc:
    :param wf_corr_num_proc:
    :param rpa_num_quad_points:
    :param size_freq_integ_group:
    :param corr_occ:
    :param corr_virt:
    :param ev_sc_iter:
    :return:
    """
    WF_CORRELATION1 = xc.WF_CORRELATION_add()
    WF_CORRELATION1.Method = 'RI_RPA_GPW'
    WF_CORRELATION1.Eri_method = 'OS'  # for non-periodic
    WF_CORRELATION1.Number_proc = wf_corr_num_proc  # -1 to use all; 16 in the ref paper
    WF_CORRELATION1.Memory = max_memory_wf
    ##### RI RPA ######
    RI_RPA = WF_CORRELATION1.RI_RPA
    RI_RPA.Rpa_num_quad_points = rpa_num_quad_points
    RI_RPA.Size_freq_integ_group = size_freq_integ_group
    RI_RPA.Gw = 'TRUE'
    ###### HF ######
    RI_RPA_HF = RI_RPA.HF_add()
    RI_RPA_HF.Fraction = 1.0
    RI_RPA_HF.SCREENING.Eps_schwarz = eps_schwarz
    RI_RPA_HF.SCREENING.Screen_on_initial_p = 'FALSE'
    RI_RPA_HF.MEMORY.Max_memory = max_memory_hf
    ###### RI_G0W0 ######
    RI_G0W0 = RI_RPA.RI_G0W0
    RI_G0W0.Corr_occ = corr_occ
    RI_G0W0.Corr_virt = corr_virt
    RI_G0W0.Ev_sc_iter = ev_sc_iter
    RI_G0W0.Analytic_continuation = 'PADE'
    #RI_G0W0.Fermi_level_offset = 0.1  #  this was a serious problem. put to default
    RI_G0W0.Crossing_search = 'NEWTON'
    RI_G0W0.Ri_sigma_x = '.TRUE.'  # x with RI: very important!
