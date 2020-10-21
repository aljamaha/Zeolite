import os, json
from ase import io
import pandas as pd

'''Prepare Input data for PCA Analyasis'''

df = pd.DataFrame()

with open('/home/aljama/BEA/Pd2/data/data.json','r') as read:
	data = json.load(read)

with open('/home/aljama/BEA/Pd2/data/data_output.json','r') as read:
	data_output = json.load(read)

def Al_Al_distance(atoms):
	'''
	Inputs: atoms object
	Output: Al-Al pair distance
	'''
	n_Al, Al_index = 0,[] #number of Al atoms, index of Al atoms
	for atom in atoms:
		if atom.symbol == 'Al':
			n_Al += 1
			Al_index.append(atom.index)
	if n_Al == 2:
		distance = atoms.get_distance(Al_index[0],Al_index[1])
	else:
		distance = 0

	return round(distance,3)

def n_atoms(atoms, atom_name):
	'finds the number of Si or O in an atom'
	n = 0
	for atom in atoms:
		if atom.symbol == atom_name:
			n+=1

	return n

def Pd_O(atoms):
	'finds the number of Si or O in an atom'
	
	distances = []
	for atom in atoms:
		if atom.symbol == 'Pd':
			Pd_index = atom.index

	for atom in atoms:
		if atom.symbol == 'O':
			distances.append(atoms.get_distance(Pd_index, atom.index))
	return sum(distances)/len(distances)
		
'''
def output_energy(item):
	'finds the energy of the item in data.json'
	for calc in data_output:
		if 'sp' in calc:
			print(calc)
			if data_output[calc]['status'] == 'complete':
				if data_output[calc]['original_info']['reference'] == item:
						print(calc, item)
						E = data_output[calc]['energy']
						break
	return E
'''

'optimum Pd2 structuresa'
Pd2 = [1486, 2078, 925, 820, 327, 547, 1778, 588, 292, 772, 1250, 314, 785, 356, 1041, 1205, 283, 706, 1034, 1023, 1702, 1417, 341, 1279, 803, 400, 463, 1002, 526, 566, 470, 604, 1862, 1525, 376, 1276, 1434, 435, 324, 874, 1533, 931, 687, 347, 956, 1819, 613, 754, 306, 413, 1053, 484, 358, 1560, 1287, 294, 765, 1914, 887, 429, 1304, 860, 723, 1548, 438, 1897, 1853, 496, 493, 943]

'calculations folders'
folders = os.system("ls /home/aljama/BEA/Pd2/calculations  > tmp")
folders = [line.rstrip('\n') for line in open('tmp')]
os.system('rm tmp')

df = {}
'initiate data frame'
df['Al-Al MR'] = []
df['Al 4MR']   = []
df['Al 5MR']   = []
df['Al 6MR']   = []
df['O']        = []
df['Si']       = []
df['Al-Al distance'] = []
df['Pd-O']       = []
df['energy']   = []

'Prepare PCA data'
PCA_data = {}
for item in Pd2:

	PCA_data[item] = {}
	
	'identify correct folder'
	for f in folders:
		if f[0:f.find('-')]+'-' == str(item)+'-' and 'sp' in f and 'omega' in f:

			atoms = io.read('/home/aljama/BEA/Pd2/calculations/'+f+'/qm-initial.traj')
			'''
			'Al-Al MR'
			if data[str(item)+'.traj']['Al-Al MR'] == []:
				PCA_data[item]['Al-Al MR'] = 0
			else:
				PCA_data[item]['Al-Al MR'] = data[str(item)+'.traj']['Al-Al MR'][0]

			'Al MR'
			PCA_data[item]['Al 4MR'] = data[str(item)+'.traj']['Al MR']['4']
			PCA_data[item]['Al 5MR'] = data[str(item)+'.traj']['Al MR']['5']
			PCA_data[item]['Al 6MR'] = data[str(item)+'.traj']['Al MR']['6']

			'Al-Al distance'
			PCA_data[item]['Al-Al distance'] = Al_Al_distance(atoms)

			'# of Si and O atoms'
			PCA_data[item]['O'] = n_atoms(atoms, 'O')
			PCA_data[item]['Si'] = n_atoms(atoms, 'Si')

			'Avg. Pd-O distance'
			PCA_data[item]['Pd-O'] = Pd_O(atoms)

			'''

			'DataFrame'
			if data[str(item)+'.traj']['Al-Al MR'] == []:
				df['Al-Al MR'].append(0)
			else:
				 df['Al-Al MR'].append(data[str(item)+'.traj']['Al-Al MR'][0])
			
			df['Al 4MR'].append(data[str(item)+'.traj']['Al MR']['4'])
			df['Al 5MR'].append(data[str(item)+'.traj']['Al MR']['5'])
			df['Al 6MR'].append(data[str(item)+'.traj']['Al MR']['6'])
			df['Al-Al distance'].append( Al_Al_distance(atoms) )
			df['O'].append( n_atoms(atoms, 'O'))
			df['Si'].append( n_atoms(atoms, 'Si'))
			df['Pd-O'].append(Pd_O(atoms))
			#df['energy'].append(output_energy(item))

#new_df = pd.DataFrame(df)
with open('PCA_data.json','w') as write_file:
    json.dump(df, write_file, indent=4)
