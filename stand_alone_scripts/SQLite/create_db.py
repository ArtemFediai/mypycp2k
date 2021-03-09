"""
create the table where to put simulation results
InchI is used as unique molecule:
https://www.rdkit.org/docs/source/rdkit.Chem.inchi.html

"""
import sqlite3
from db.db import Molecule

def main():
    conn = sqlite3.connect("molecule.db")
    curs = conn.cursor()

    curs.execute("""CREATE TABLE IF NOT EXISTS simulations (
                mol_id integer,
                InChI text,
                mol_name text,
                SMILE text,
                num_atoms integer,
                num_orb integer,
                basis_set text,
                date text,
                homo real,
                lumo real,
                gw_homo real,
                gw_lumo real
                );""")

    # example, how to insert a mol into db table 'simulations'
    my_molecule = Molecule(mol_name='hexomethane')
    my_molecule.insert_to_db(conn, 'simulations')
    conn.commit()
    conn.close()

    print("I am done!")

if __name__ == '__main__':
    main()
