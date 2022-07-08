def generate_esp_point_charges(FORCE_EVAL,
                               r_min_scale=1.4,
                               r_max_scale=2.0,
                               print_esp_charges=True,
                               filename=None,
                               vdw_radii_table='CAMBRIDGE'):
    """
    fit point charges. by default between 1.4 and 2 vdW radii
    @param vdw_radii_table: type of the mentioned Table
    @param print_esp_charges: wheather to print esp
    @param FORCE_EVAL: pyscf object
    @param r_max_scale: rmax = R_MIN_SCALE * vdW_radius
    @param r_min_scale: rmin = R_MAX_SCALE * vdW_radius
    @return: None
    """

    PROPERTIES = FORCE_EVAL.PROPERTIES
    RESP = PROPERTIES.RESP

    RESP.Restrain_heavies_to_zero = '.FALSE.'
    # RESP.Use_repeat_method = use_repeat_method   #<-- only for periodic:(

    SPHERE_SAMPLING = RESP.SPHERE_SAMPLING
    SPHERE_SAMPLING.Auto_vdw_radii_table = vdw_radii_table
    SPHERE_SAMPLING.Auto_rmin_scale = r_min_scale
    SPHERE_SAMPLING.Auto_rmax_scale = r_max_scale

    if print_esp_charges:
        PRINT = RESP.PRINT
        RESP_CHARGES_TO_FILE = PRINT.RESP_CHARGES_TO_FILE
        RESP_CHARGES_TO_FILE.Section_parameters = 'LOW'

