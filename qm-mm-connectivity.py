#import packages
import sys
from molmod import *

'input format: python qm-mm-connectivity'

'raise an error when less than two arguments are provided'
if len(sys.argv) !=3:
  sys.stderr.write(sys.argv[0]+"\n")
  sys.stderr.write("*************************************************************************************************** \n")
  sys.stderr.write("Add forcefield atom types and connectivity to XYZ coordinates for zeolite QM/MM calculations with Q-Chem.\n\n")
  sys.stderr.write("Usage: \n")
  sys.stderr.write("qmmm-connectivity-standalone.py <XYZ file> <number of QM atoms (must be at the top of the XYZ file)>\n")
  sys.stderr.write("*************************************************************************************************** \n")
  quit()

mol = Molecule.from_file(sys.argv[1])

mol.set_default_masses()

assert(mol.graph is None)

mol.set_default_graph()

qm_atoms = range(0,int(sys.argv[2]))

link_O = []

print(mol.graph.neighbors)

for i, ns in mol.graph.neighbors.iteritems():
    if i in qm_atoms:
        indexes = [n for n in ns]
        for j in range(0,len(indexes)):
                if not indexes[j] in qm_atoms:
                   link_O.append(indexes[j])

link_Si = []
for i, ns in mol.graph.neighbors.iteritems():
    if i in link_O:
        indexes = [n for n in ns]
        for j in range(0,len(indexes)):
                if not indexes[j] in qm_atoms and mol.symbols[indexes[j]] == 'Si':
                   link_Si.append(indexes[j])

for i, ns in mol.graph.neighbors.iteritems(): 
    indexes = [n for n in ns]
    if i in qm_atoms:
         type = 0
    elif i in link_Si:
         type = -24
         for j in range(0,len(indexes)):
                if indexes[j] in link_O:
                          type += 1
    elif mol.symbols[i] == 'O':
        type = -1
    elif mol.symbols[i] == 'Si' or mol.symbols[i] == 'Al':
        type = -2
    elif mol.symbols[i] == 'H':
        type = -3
    print("%s\t % f\t % f\t % f\t % i\t %i\t %i\t %i\t %i\t" % (str(mol.symbols[i]), mol.coordinates[i][0]/angstrom, mol.coordinates[i][1]/angstrom, mol.coordinates[i][2]/angstrom, type, indexes[0]+1 if len(indexes) > 0 else 0, indexes[1]+1 if len(indexes) > 1 else 0, indexes[2]+1 if len(indexes) > 2 else 0, indexes[3]+1 if len(indexes) > 3 else 0))


