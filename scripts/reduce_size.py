import os, json
from ase import io

'''Reduce size of opt.out files for storage purposes'''

'Inputs'
n_lines   = 76000	#number of preserved lines at the beginning and end of opt.out'
n_minimum = 300000	#minimum number of lines to be eligible for reduction
calc_dir  = '/home/aljama/BEA/H/calculations'

def reduce_size(n_line, n_minimum):
	'rewrites opt.out by stripping many lines and maintaining important ones'
	'n_lines: number of preserved lines at the beginning and end of opt.out'

	'opening files and over-writing'
	f = open('opt.out')
	lines = f.read().split("\n")
	f.close()

	'inputs'
	low  = n_lines
	L =  len(lines)
	high = L-low
	
	if L > n_minimum:
		'writing new file'
		f = open('opt2.out','w')
		for index, line in enumerate(lines):
			try:
				if index < low:
					f.write(line+'\n')
				elif index > high:
					f.write(line+'\n')
				elif '  Gradient  ' in line or '  Displacement ' in line or "Energy change" in line:
					f.write(line+'\n')
				elif 'Total energy ' in line:
					f.write(line+'\n')
			except:
				pass

	f.write('original opt.out file is trimmed to reduce file size')
	f.close()

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

def print_qm(wd):
	'print traj file of qm atoms'
	try:
		n_Al,n_Si, n_O = Al_Si_atoms()	#number of Al/Si atoms in qm region
		H_qm = H_qm_region(n_O, n_Si, n_Al)  #H added to qm region to justify MM region
		data = data_load()
		n_qm = len(data['qm_region']) + H_qm #total atoms in qm region
		n_qm = str(n_qm)
		print(n_qm)
		os.system('cp ~/scripts convert-qchem-output-to-ase.py .')
		os.system('python convert-qchem-output-to-ase.py '+n_qm)
		if '1.xyz' in folders_under(wd+'/traj') == True:
			tmp  = True
		else:			
			tmp = False
	except:
		tmp = False

	return tmp

def data_load():
	'Load dir_data.json'
	with open('dir_data.json','r') as read:
		data = json.load(read)

	return data

def qm_traj_print(wd):
	'saves traj files of the qm region'

	os.chdir(wd)
	'Load dir_data.json'
	try:
		with open('dir_data.json','r') as read:
			data = json.load(read)
			qm   = data['QMMM length']
		os.system('cp ~/scripts/convert-qchem-output-to-ase.py .')
		os.system('python convert-qchem-output-to-ase.py '+str(qm))
		if '1.xyz' in folders_under(wd+'/traj') == True:
			tmp = 'True'
		else:
			tmp = 'False'
			print('WARNING! traj folder could not be created')
	except:
		tmp = 'False'
		print('WARNING! traj folder could not be created')

	return tmp

def folders_under(wd):
	'returns files and folders in the wd'
	os.chdir(wd)
	os.system("ls > tmp")
	folders = [line.rstrip('\n') for line in open('tmp')]  #this is the best way
	return folders

'Start here ..'
os.chdir(calc_dir)
folders = folders_under(calc_dir)
folders = ['1216-opt-omegab97x-d-def2-svp-ref-61.traj']

for folder in folders:
	if 'opt' in folder:
		print(folder)
		os.chdir(calc_dir+'/'+folder)
		tmp = print_qm(calc_dir+'/'+folder)
		print(tmp)
		if tmp == True:
			try:
				reduce_size(n_lines, n_minimum)
			except:
				print('could not reduce size!')
		else:
				print('could not reduce size!')
