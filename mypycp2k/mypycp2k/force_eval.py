def generate_esp_point_charges(FORCE_EVAL,
                               r_min_scale=1.4,
                               r_max_scale=2.0,
                               print_esp_charges=True,
                               filename=None,
                               vdw_radii_table='CAMBRIDGE',
                               stride=1,
                               elements=None,
                               use_tm_radii=False,
                               integer_total_charge=True,
                               restrain_heavies_to_zero=False
                               ):
    """
    fit point charges. by default between 1.4 and 2 vdW radii
    @param use_tm_radii: use Turbomole vdW radii defined below
    @param elements: list of elements from xyz file
    @param vdw_radii_table: type of the mentioned Table
    @param print_esp_charges: wheather to print esp
    @param FORCE_EVAL: pyscf object
    @param r_max_scale: rmax = R_MIN_SCALE * vdW_radius
    @param r_min_scale: rmin = R_MAX_SCALE * vdW_radius
    @return: None
    """
    
    # TM dict. vdW radii
    dict_au = {
        'H': 2.267671,
        'C': 3.212534,
        'N': 2.929075,
        'O': 2.872383,
        'F': 2.777897,
        'Al': 3.590479,
    }

    from scipy.constants import physical_constants
    A_in_au = physical_constants['Bohr radius'][0]*1E10
    dict_A = {key: value*A_in_au for (key, value) in dict_au.items()}  # in A

    PROPERTIES = FORCE_EVAL.PROPERTIES
    RESP = PROPERTIES.RESP

    RESP.Integer_total_charge = integer_total_charge
    RESP.Restrain_heavies_to_zero = restrain_heavies_to_zero
    # RESP.Use_repeat_method = use_repeat_method   #<-- only for periodic:(
    RESP.Stride = stride
    SPHERE_SAMPLING = RESP.SPHERE_SAMPLING

    SPHERE_SAMPLING.Rmin_kind = []
    SPHERE_SAMPLING.Rmax_kind = []
    if use_tm_radii:
        for element in elements:
            SPHERE_SAMPLING.Rmin_kind.append(' '.join([str(r_min_scale*dict_A[element]), element]))
            SPHERE_SAMPLING.Rmax_kind.append(' '.join([str(r_max_scale*dict_A[element]), element]))

    else:
        SPHERE_SAMPLING.Auto_vdw_radii_table = vdw_radii_table
        SPHERE_SAMPLING.Auto_rmin_scale = r_min_scale
        SPHERE_SAMPLING.Auto_rmax_scale = r_max_scale

    if print_esp_charges:
        PRINT = RESP.PRINT
        RESP_CHARGES_TO_FILE = PRINT.RESP_CHARGES_TO_FILE
        RESP_CHARGES_TO_FILE.Section_parameters = 'LOW'

