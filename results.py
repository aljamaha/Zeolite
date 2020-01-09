import json, os
import matplotlib.pyplot as plt
from ase import atoms, io

'''
analyzes the results of completed calcluations
'''

'Inputs'
plotting    	    = False		#if True, plot results for each reference structure
plotting_overall    = False		#if True, make an overall plot of all results
data_dir    = '/home/aljama/CHA/data/'			#dir where json data are saved
calc_dir    = '/home/aljama/CHA/calculations/'		#dir where calculations are done
results_dir = '/home/aljama/CHA/results-analysis/' 	#dire where results are to be saved

'Load data from json files'
with open(data_dir+"data_output.json", "r") as read_file:
    data_output = json.load(read_file)

with open(data_dir+"data.json", "r") as read_file:
    data_original = json.load(read_file)

def calc_index(index):
	'''
	finds the entry in data_output and matches it to the item in data_original ['index']
	Input : index in data_original
	Output: entry in data_output
	'''
	output = 'none' #make default is none, uncless calculation is available

	for item in data_output:
		if data_output[item]['index'] == index:
			stat = '-sp-' in item
			if stat == True:
				output = item
	return output

def Al_Al(atoms):
	'''
	finds distance between 2 Al atoms in qm region
	Input : ase atoms object
	Output: Al-Al distance [unless there is one Al, distance then is 0] and number of Al atoms
	'''

	n_Al = 0 		#number of Al atoms
	Al_index = []	#index of Al atom

	for atom in atoms:
		if atom.symbol == 'Al':
			n_Al += 1
			Al_index.append(atom.index)

	if n_Al == 1:
		distance = 0
	elif n_Al == 2:
		distance = atoms.get_distance(Al_index[0],Al_index[1])
	
	return distance, n_Al

def O_O(atoms, n_Al):
	'''
	finds the distance between two oxygen atoms where H is bonded
	Inputs:
		atoms: ase atoms object
		n_Al : number of Al atoms in qm region
	Output: distance between two oxygen atoms where H is bonded [0 if n_Al = 1]
	'''

	H_index, O_index = [],[]  #indexes of H and O atoms

	if n_Al == 1:
		distance = 0
	elif n_Al == 2:
		for atom in atoms:
			'find H atoms in qm region'
			if atom.symbol == 'H':
				H_index.append(atom.index)

		for H in H_index:
			'find closest O to the H'
			d_min = 100 #start with a large value (dummy value)

			for atom in atoms:
				if atom.symbol == 'O':
					if atoms.get_distance(H, atom.index) < d_min:
						d_min  =  atoms.get_distance(H, atom.index)
						O_atom = atom.index

			O_index.append(O_atom)
	
		distance = atoms.get_distance(O_index[0], O_index[1])
	
	return distance

'accumulate reference entries (templates from which calculations were created and run)'
references = {} #references for data

for item in data_original:
	'entries in references dictionary'
	if data_original[item]['reference'] not in references:
		references[data_original[item]['reference']] = []

for ref in references:
	'name of folders that share same reference'
	for item in data_original:
		if data_original[item]['reference'] == ref:
			if ref != item:
				references[ref].append(item)

'accumulate traj files'
O_d_all, E_all = [],[] #those are for all calculations to get a comprehensive picture

for ref in references:
	'loop over each reference'
	x_pos, E, first_item, O_d, label = [],[], True, [],[]	#For plotting purposes

	if os.path.isdir(results_dir+ref) == False:
		'create folder if it does not exist'
		os.system('mkdir '+results_dir+ref)

	os.chdir(results_dir+ref)

	for item in references[ref]:
		'each item under reference'
		index = item[0:-5] #remves .traj from the name
		data_output_entry = calc_index(index) #check corresponding name in data_output
		if data_output_entry != 'none':
			'check calcuation dir is available'
			if data_output[data_output_entry]['status'] == 'complete':
				'check calc is completed, then copy traj files to new folder'
				os.system('cp '+calc_dir+'/'+data_output_entry+'/qm-initial.traj '+item[0:-5]+'.traj')
				x_pos.append(int(index))	#x-asis position
				if first_item == True:
					E_ref = data_output[data_output_entry]['energy']
				E.append( (data_output[data_output_entry]['energy']- E_ref)*27.2114 ) #convert from Hartree to e. values of energies for y-axis
				E_all.append( (data_output[data_output_entry]['energy']- E_ref)*27.2114 ) 
				first_item = False

				'Al-Al distance'
				atoms = io.read(calc_dir+data_output_entry+'/qm-initial.traj')
				Al_distance, n_Al = Al_Al(atoms)
				data_output[data_output_entry]['Al-Al distance'] = round(Al_distance,3)

				'O-O distance'
				O_O_distance = O_O(atoms, n_Al)
				data_output[data_output_entry]['O-O distance'] = round(O_O_distance,3)
				O_d.append( round(O_O_distance,3) )
				O_d_all.append( round(O_O_distance,3) )
				
				label.append(index)
				
	if plotting == True:
		plt.bar(x_pos, E, align='center', alpha=1)
		plt.xticks(x_pos, x_pos)
		plt.ylabel('Energy (eV)')

		plt.show()
	
		plt.plot(O_d, E, 'o', markersize=6)
		plt.xlabel('O-O Distance (A)', fontsize = 10)
		plt.ylabel('Energy (eV)', fontsize = 10)
		for i, lab in enumerate(label):
			plt.text(O_d[i], E[i], lab)
		plt.show()

if plotting_overall == True:
	plt.plot(O_d_all, E_all, 'o', markersize=6)
	plt.xlabel('O-O Distance (A)', fontsize = 10)
	plt.ylabel('Energy (eV)', fontsize = 10)
	plt.xlim([4,12])
	plt.show()

'save data with Al-Al distance'
with open(data_dir+"data_output.json", "w") as write_file:
    json.dump(data_output, write_file, indent=4)
