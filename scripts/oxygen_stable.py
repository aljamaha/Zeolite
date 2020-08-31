from ase import io
from molmod import *

atoms = io.read('qm-initial.traj')
atoms.write('tmp.xyz')

def individual_NL(index, N_list):
	'creates a neighboring list for atom of interest (needs only atom index and complete NL'
	n_list = {} #local neighbor list specific to Si atoms next to O next to terminal Si
	n_list['Si'], n_list['O'] = {},{}
	n_list['Si']['N'], n_list['O']['N'] = [],[]
	n_list = identify_N(N_list[index], N_list, n_list, index)

	return n_list

def neighbor_list(xyz_file):
	'''
	Inputs:  xyz coordinates
	Outputs: dictionary of the neighbor list (for all atoms)
		e.g. [atom index]: [neighboring atoms]
	'''
	mol = Molecule.from_file(xyz_file)
	mol.set_default_masses()
	mol.set_default_graph()
	return mol.graph.neighbors

def identify_N(N_list_Al, N_list, neighbors, Al):
	'''
	Identify Al neighboring Si and O
	Inputs :
		N_list[Al]: neighbor list of Al
		N_list    : dictionary of the neighboring list of all atoms
		neighbors : dictionary of neighbors (Si/O)[N,NN,NNN]
		Al 	  : element number of first Al
	Outputs:
		neighbors: Si and O neighboring Al
	'''
	for item in N_list_Al:
		'loop over O attached to Si'
		neighbors['O']['N'].append(item)
		'''
		for neighbor in N_list[item]:
			'loop over neighbors of Si attached to neighbor'
			if neighbor == Al:
				continue
			elif neighbor in neighbors['Si']['N']:
				continue
			else:
				neighbors['Si']['N'].append(neighbor)
		'''
	return neighbors

N_list = neighbor_list('tmp.xyz')
oxygen_N_Al    = []

for atom in atoms:
	if atom.symbol == 'Al':
		tmp = individual_NL(atom.index, N_list)
		print(tmp)
		
		
