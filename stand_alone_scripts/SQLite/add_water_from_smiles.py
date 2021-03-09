"""
Suppose we have something saved as SMILEs (in this example: water molecule)
This scripts adds it into the table
"""

from rdkit import Chem

m = Chem.MolFromSmiles('O')

InChI = Chem.MolToInchi(m)
InChIKey = Chem.MolToInchiKey(m)


print("Inchi: ", InChI)
print("Inchi key: ", InChIKey)
