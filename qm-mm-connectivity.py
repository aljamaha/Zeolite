#import packages
import sys
from molmod import *

'''
Objective: convert an xyz coordinates to qmmm molecule format (usign neighbors)
Input format: python qm-mm-connectivity <xyz coordinates> <# of qm atoms>
Output: $molecule section in qm-mm
'''

if len(sys.argv) !=3:
	'raise an error when less than two arguments are provided'
	sys.stderr.write("Need two arguments in the input [<xyz format> <#of qm atoms>]\n")
	quit()

mol = Molecule.from_file(sys.argv[1])	#reads the input <xyz coordinates>
mol.set_default_masses()		#reads the mass of each atom in xyz input
assert(mol.graph is None)
mol.set_default_graph()
qm_atoms = range(0,int(sys.argv[2]))	#number of qm atoms based on input file

link_O = []				#O-atoms connected to qm atoms
#for i, ns in mol.graph.neighbors.iteritems():
for 
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


