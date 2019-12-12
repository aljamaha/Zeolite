import os, json

'Analyze the output of qmmm calculations and exports data to output_data.json'

'Inputs'
calc_dir = '/home/aljama/CHA/calculations/'	#directory where caluculatiosn are saved
data_dir = '/home/aljama/CHA/data/'		#directory where data are saved
scripts_dir = '/home/aljama/scripts/'

cwd = os.getcwd()

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

	step, energy,steps = 0,[],[]

	for line in lines:
		'read energy per step'
		if 'Total energy in the final basis set' in line:
			step += 1
			if step == 1:
				E0 = float(line[-14:])
			energy.append(float(line[-14:])- E0)
			steps.append(step)

	Final_Energy = energy[-1] + E0

	return Final_Energy

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

folders = folders_list(calc_dir)	#folders in calculations/ directory

'open saved json file'
with open(data_dir+"/data_output.json", "r") as read_file:
	data = json.load(read_file)		#output data

with open(data_dir+"/data.json", "r") as read_file:
	data_original = json.load(read_file)	#original data details

'data extraction'
for folder in folders:
	tmp_status = '' #temporary set that calc is not complete
	if folder in data:
		if 'status' in data[folder]:
			if data[folder]['status'] == 'complete':
				'if it passes those tests, calc is complete and no need for further check'
				tmp_status = 'complete'
	if tmp_status == 'complete':
		'completed calc has been stored'
		continue
	else:
		'if calculation is incomplete, then proceed to check status'
		data[folder] = {}
		os.chdir(calc_dir+'/'+folder)

		for index, item  in enumerate(folder):
			if item == '-':
				'extracts index (used to extract original data)'
				data[folder]['index'] = folder[0:index]
				ref = str(folder[0:index])+'.traj'	#reference in original calculation
				break

		'extract calculation detail'
		data[folder]['calc']   = folder[index+1:]
		data[folder]['status'] = calc_status(os.getcwd())[0]

		if data[folder]['status'] == 'complete':
			'extract energy if calculation is complete'
			output_file = calc_status(os.getcwd())[1]
			E = extract_energy(output_file)
			data[folder]['energy'] = E

			'print traj of the qm region'
			os.system('python '+scripts_dir+'/convert-qchem-output-to-ase.py '+output_file+' '+'input.xyz'+' '+str(len(data_original[ref]['qm_region'])+16))

			'print full traj of all atoms'
			os.system('python '+scripts_dir+'qchem-to-ase-all-atoms.py '+output_file+' '+'input.xyz'+' '+str(data_original[ref]['total_atoms']))

for item in data:
	if data[item]['status'] == 'incomplete':
		print(item, data[item])

'saving output data'
with open(data_dir+"/data_output.json", "w") as write_file:
    json.dump(data, write_file, indent=4)
