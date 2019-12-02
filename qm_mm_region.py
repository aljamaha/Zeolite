#!/Users/hassanaljama/opt/anaconda3/bin/python

from ase import io
import pickle
from functions import identify_N


def qm_mm_region(data,traj, struc_dir, cwd, N_list):
	'''
	identifies elements in qm region of a zeolite structure based on Al atoms
	inputs : 
		 data  : data dictionary
		 traj  : name of the traj of the structure
		 cwd   : current working directory
		 N_list: neighbor list of all atoms
	outputs: data dictionary
	'''

	data[traj]['qm_region'] =  []			#add entry for qm_region
	atoms     = io.read(cwd+'/'+struc_dir+'/'+traj)
	neighbors, neighbors['Si'], neighbors['O'] = {},{'N':[]},{'N':[]}

	'identify Al atoms'
	Al_atoms = [] 	
	for index, atom in enumerate(atoms):
		if atom.symbol == 'Al':
			Al_atoms.append(index)

	'identify neighboring Si/O'
	for Al in Al_atoms:
		data[traj]['qm_region'].append(Al)
		neighbors = identify_N(N_list[Al], N_list, neighbors, Al)
	
	'add Si/O elements to qm region'
	for item in neighbors:
		if neighbors[item]['N'] in data[traj]['qm_region']:
			continue
		else:
			data[traj]['qm_region'].append(neighbors[item]['N'])

	return data
	
data = {}
N_list    = pickle.load( open( "save.p", "rb" ) )
qm_mm_region(data,traj, struc_dir, cwd, N_list):

