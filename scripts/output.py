import os, json, sys
from ase import io
from copy import deepcopy

'Analyze the output of qmmm calculations and exports data to output_data.json'

'Inputs'
dir_name     = 'BEA/Pd1'
traj         = False   #prints traj files of qm region
full	     = False  #prints traj files of full atoms

'Directroies'
calc_dir    = '/home/aljama/'+dir_name+'/calculations/'	#directory where caluculatiosn are saved
data_dir    = '/home/aljama/'+dir_name+'/data/'		#directory where data are saved
scripts_dir = '/home/aljama/scripts/'		#directory where zeolites scripts are

def folders_list(wd):
	'list of folders in a directory'
	os.chdir(wd)
	os.system("ls > tmp")
	folders = [line.rstrip('\n') for line in open('tmp')]
	folders.remove('tmp')
	os.system('rm tmp')

	return folders

def extract_energy(output_file):
	'''
	Extract output energy
	outputs: Energy of the structure in the output file
	'''
	f = open(output_file,'r')
	lines = f.read().split("\n")
	f.close()

	E = 'nan'

	for line in lines:
		'read energy per step'
		if 'Total energy' in line:
			ind = line.find('=')
			E   = float(line[ind+1:])

	return E

def calc_status(wd):
	'''
	check calcualtion status
	Output: 
		status of calculation (complete/incomplete)
		name of output folder
	'''
	folders = folders_list(wd)
	output_folder,status = 'None','incomplete'	#initially assign none until one is found

	for folder in folders:
		if "out" in folder:
			output_folder = folder
			f = open(folder)
			text = f.read()
			f.close()
			if 'Have a nice day' in text:
				status = 'complete'

	return status, output_folder

def Al_Si_atoms():
	'finds number of Al atoms in qm region'
	atoms = io.read('qm-initial.traj')
	n_Al, n_Si, n_O = 0, 0, 0
	for atom in atoms:
		if atom.symbol == 'Al':
			n_Al += 1
		elif atom.symbol == 'Si':
			n_Si += 1
		elif atom.symbol == 'O':
			n_O  += 1

	return n_Al, n_Si, n_O

def H_qm_region(n_O, n_Si, n_Al):
	'''identifies the number of H atoms in qm region added to justify dangling bonds
	Inputs: 
		n_O:  number of O  in qm region
		n_Si: number of Si in qm region
		n_Al: number of Al in qm region
	Outputs:
		number of additional H in qm region to prevent dangling Si bonds
	'''

	num_oxygen_bonds = n_O*2	#total number of bonds oxygen can make in QM region
	num_oxygen_bonds_to_Al = n_Al*4
	num_oxygen_bonds_to_Si = num_oxygen_bonds - num_oxygen_bonds_to_Al
	tot_num_Si_bonds_needed = n_Si*4
	tot_Si_bonds_with_H_needed = tot_num_Si_bonds_needed - 	num_oxygen_bonds_to_Si

	return tot_Si_bonds_with_H_needed 
	
folders = folders_list(calc_dir)	#folders in calculations/directory

'open saved json file'

if os.path.isfile(data_dir+"/data_output.json") == True:
	with open(data_dir+"/data_output.json", "r") as read_file:
		data = json.load(read_file)		#output data
else:
	data = {}
	with open(data_dir+"/data_output.json", "w") as write_file:
  		json.dump(data, write_file, indent=4)

with open(data_dir+"/data.json", "r") as read_file:
	data_original = json.load(read_file)	#original data details

'delete items under calculations that are not opt/sp'
folders_copy = deepcopy(folders)
for folder in folders_copy:
	if 'opt' not in folder :
		if 'sp' not in folder:
			folders.remove(folder)

'data extraction'
for folder in folders:

	print('calculation: ', folder)
	tmp_energy = '' #assign an empty initial value
	
	'check status of calculation and energy'
	if folder in data:
		if 'energy' in data[folder]:
			tmp_energy = 'extracted'
	else:
		'create new entry if it does not exist'
		data[folder] = {}	#create new entry in dictionary

	os.chdir(calc_dir+'/'+folder)

	for index, item  in enumerate(folder):
		if item == '-':
			'extracts index (used to extract original data)'
			data[folder]['index'] = folder[0:index]					
			ref = str(folder[0:index])+'.traj'	#reference in original calculation
			break

	'name of output file'
	output_file = calc_status(calc_dir+'/'+folder)[1]

	#if tmp_energy != 'extracted':
	if tmp_energy != '??':
		'extract calculation detail/Energy'
		data[folder]['calc']   = folder[index+1:] 
		data[folder]['status'] = calc_status(calc_dir+'/'+folder)[0]
		
		if data[folder]['status'] == 'complete':
			'extract energy if calculation is complete'
			E = extract_energy('opt.out')
			data[folder]['energy'] = E
			if E == 'nan':
				data[folder]['status'] = 'incomplete'
				tmp_energy = ''
			else:
				tmp_energy = 'extracted'

	if tmp_energy in ['extracted','']:	
			'print qm traj files if they do not exist'
			os.chdir(calc_dir+'/'+folder)
			os.system('cp '+scripts_dir+'/convert-qchem-output-to-ase.py .')

			if traj == True:
				if 'opt' in folder:
					n_Al,n_Si, n_O = Al_Si_atoms()	#number of Al/Si atoms in qm region
					H_qm = H_qm_region(n_O, n_Si, n_Al)  #H added to qm region to justify MM region
					n_qm = len(data_original[ref]['qm_region']) + H_qm #total atoms in qm region
					n_qm = str(n_qm)
					os.system('python convert-qchem-output-to-ase.py '+n_qm)

			'print full traj of all atoms if they do not exist'
			if full == True:
				if 'opt' in folder:
					os.system('python '+scripts_dir+'qchem-to-ase-all-atoms.py '+output_file+' '+'input.xyz'+' '+str(data_original[ref]['total_atoms']))
			
			#'print traj files of surroundings'
			#if surroundings == True:
			#	os.system('python '+scripts_dir+'qchem-to-ase-surroundings.py')

			'adding information from original calculation'
			data[folder]['original_info'] = data_original[ref]


print('Incomplete calculations:')
for item in data:
	if data[item]['status'] == 'incomplete':
		print(item)

'saving output data'
with open(data_dir+"/data_output.json", "w") as write_file:
    json.dump(data, write_file, indent=4)
