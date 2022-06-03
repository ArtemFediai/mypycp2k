def add_elements(subsys, elements, basis, aux_basis, pot):
    """
    Add chemichal elements sections
    :param subsys: pycp2k object
    :param elements: list of checmichal elements
    :param basis: basis set name
    :param aux_basis: ri axu basis set name
    :param pot: pseudopotential name
    :return: add kinds to KIND object
    """
    if isinstance(elements, list):
        elements = list(elements)

    for element in elements:
        KIND = subsys.KIND_add(element)
        KIND.Element = element
        KIND.Basis_set = basis
        KIND.Potential = pot
        KIND.Ri_aux_basis_set = aux_basis


def set_unperiodic_cell(subsys, abc):
    CELL = subsys.CELL
    CELL.Abc = abc
    CELL.Periodic = 'NONE'

def set_periodic_cell(subsys, abc):
    CELL = subsys.CELL
    CELL.Abc = abc
    CELL.Periodic = 'XYZ'

def set_topology(subsys, xyz_file_name):
    TOPOLOGY = subsys.TOPOLOGY
    TOPOLOGY.Coord_file_format = 'xyz'
    TOPOLOGY.Coord_file_name = xyz_file_name


def center_coordinates(subsys):
    TOPOLOGY = subsys.TOPOLOGY
    CENTER_COORDINATES = TOPOLOGY.CENTER_COORDINATES
    CENTER_COORDINATES.Section_parameters = ''
