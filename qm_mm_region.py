from ase import io
import pickle

'inputs'
atoms    = io.read('/Users/hassanaljama/Desktop/CHA/structures/1.traj')
N_list   = pickle.load( open( "save.p", "rb" ) )

'identify Al atoms'
Al_atoms = [] 	
for index, atom in enumerate(atoms):
	if atom.symbol == 'Al':
		Al_atoms.append(index)


