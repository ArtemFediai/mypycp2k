def generate_esp_point_charges(FORCE_EVAL,
                               r_min_scale=1.4,
                               r_max_scale=2.0,
                               print_esp_charges=True,
                               filename=None):
    """
    fit point charges. by default between 1.4 and 2 vdW radii
    @param print_esp_charges: wheather to print esp
    @param FORCE_EVAL: pyscf object
    @param r_max_scale: rmax = R_MIN_SCALE * vdW_radius
    @param r_min_scale: rmin = R_MAX_SCALE * vdW_radius
    @return: None
    """
    print('Implement here')
    PROPERTIES = FORCE_EVAL.PROPERTIES
    RESP = PROPERTIES.RESP
    SPHERE_SAMPLING = RESP.SPHERE_SAMPLING
    SPHERE_SAMPLING.Auto_vdw_radii_table = 'CAMBRIDGE'
    SPHERE_SAMPLING.Auto_rmin_scale = r_min_scale
    SPHERE_SAMPLING.Auto_rmin_scale = r_max_scale

    if print_esp_charges:
        PRINT = RESP.PRINT
        RESP_CHARGES_TO_FILE = PRINT.RESP_CHARGES_TO_FILE
        RESP_CHARGES_TO_FILE.Section_parameters = 'LOW'

