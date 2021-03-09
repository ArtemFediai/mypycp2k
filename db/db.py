class Molecule:
    """
    object, whose attrs are a coulomn in simulation table of molecule.db
    """
    def __init__(self,
                 mol_id=None,
                 InChI=None,
                 mol_name=None,
                 SMILE=None,
                 num_atoms=None,
                 num_orb=None,
                 basis_set=None,
                 date=None,
                 homo=None,
                 lumo=None,
                 gw_homo=None,
                 gw_lumo=None):
        self.mol_id = mol_id
        self.InChI = InChI
        self.mol_name = mol_name
        self.SMILE = SMILE
        self.num_atoms = num_atoms
        self.num_orb = num_orb
        self.basis_set = basis_set
        self.date = date
        self.homo = homo
        self.lumo = lumo
        self.gw_homo = gw_homo
        self.gw_lumo = gw_lumo

    def insert_to_db(self, cunnection_, table_name_):
        curs = cunnection_.cursor()
        my_tuple = self.as_tuple()
        num_col_in_db = len(my_tuple)
        test_line = f"INSERT INTO {table_name_} VALUES ({'? , '*(num_col_in_db - 1) + ' ?'})"  # last '?' w/o comma
        curs.execute(test_line, my_tuple)


    def as_tuple(self):
        dict_all_attr = self.__dict__
        return tuple(dict_all_attr[key_] for key_ in dict_all_attr if not key_.startswith('_'))
