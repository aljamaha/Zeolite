#!/Users/hassanaljama/opt/anaconda3/bin/python

import os
from ase import io
import json

'''
Running calculations on selected strucutres. Provide Inputs below
'''
'Inputs'
dir_name     = 'BEA/H'	#name of the parent dir
multiplicity = 1		#multiplicity of the structure
exchange     = 'omegab97x-d'		#'omegab97x-d' or 'B97-D3'
job_type     = 'sp' 		#either sp or opt
zeolite      = 'BEA'		#zeolite name
#calc         = [1798,1798+32]		#if continous list, input first and last. Otherwise, input individual entries
#calc	     = [506, 507, 508, 509, 510, 511, 512, 513, 514, 515, 516, 517, 518, 519, 520, 521, 522, 523, 524, 525, 526, 527, 528, 529, 530, 531, 532, 533, 534, 535, 536, 537, 3302, 3303, 3304, 3305, 3306, 3307, 3308, 3309, 3310, 3311, 3312, 3313, 3314, 3315, 3316, 3317, 3318, 3319, 3320, 3321, 3322, 3323, 3324, 3325, 3326, 3327, 3328, 3329, 3330, 3331, 3332, 3333, 4906, 4907, 4908, 4909, 4910, 4911, 4912, 4913, 4914, 4915, 4916, 4917, 4918, 4919, 4920, 4921, 4922, 4923, 4924, 4925, 4926, 4927, 4928, 4929, 4930, 4931, 4932, 4933, 4934, 4935, 4936, 4937]
#calc = [394, 395, 396, 397, 398, 399, 400, 401, 402, 403, 404, 405, 406, 407, 408, 409, 1798, 1799, 1800, 1801, 1802, 1803, 1804, 1805, 1806, 1807, 1808, 1809, 1810, 1811, 1812, 1813, 2602, 2603, 2604, 2605, 2606, 2607, 2608, 2609, 2610, 2611, 2612, 2613, 2614, 2615, 2616, 2617]
calc = [278, 279, 280, 281, 746, 747, 748, 749, 1214, 1215, 1216, 1217, 1714, 1715, 1716, 1717, 2214, 2215, 2216, 2217, 2682, 2683, 2684, 2685, 3150, 3151, 3152, 3153, 3634, 3635, 3636, 3637, 4134, 4135, 4136, 4137]
print('Directory name', dir_name)

#cont_check = input('Press Enter to Continue ..')

'Inputs (rarely need a change)'
struc_dir        = '/home/aljama/'+dir_name+'/structures_saved'
data_dir         = '/home/aljama/'+dir_name+'/data'
scripts_dir      = '/home/aljama/scripts/' 	#this is where qm-initial structure script is
calc_dir         = '/home/aljama/'+dir_name+'/calculations/'
create_input_dir = '/home/aljama/BEA/original/create_input/'

'basis set depending on job_type'
if job_type == 'sp':
	basis    	= 'def2-tzvpd' 	#basis set [def2-sv(p) or def2-tzvpd]
elif job_type == 'opt':
	basis    	= 'def2-sv(p)' 	#basis set [def2-sv(p) or def2-tzvpd]
details  = job_type+'-'+exchange+'-'+basis.replace('(','').replace(')','') #naming dir (uniqueness)

'Calculations list'
calculations = []
if len(calc) == 2:
	'continous list from start to end'
	for j in range(calc[0],calc[1]):
		calculations.append(str(j))
else:
	'individual list'
	for j in calc:
		calculations.append(str(j))

'Load data'
with open(data_dir+"/data.json", "r") as read_file:
    data = json.load(read_file)
with open(data_dir+"/data_output.json", "r") as read_file:
    data_output = json.load(read_file)

def rem_section():
	'writes details of rm section'
	g = open(create_input_dir+'/text-rm.txt','r')
	text_rm = g.read()
	g.close()
	f.write('$rem   \n')
	f.write('jobtype \t'+job_type+'\n')
	f.write('exchange   '+exchange+'\n')
	f.write('basis   \t'+basis+'\n')
	f.write('AIMD_FIXED_ATOMS \t'+str(len(fixed_atoms))+'\n')
	f.write(text_rm)
	f.write('model_system_mult '+str(multiplicity)+'\n')
	f.write('$end\n\n')

def qm_atoms_section(qm_atoms):
	'writes details of qm_atoms section'
	f.write('$qm_atoms\n')
	qm_atoms.sort() #if unsorted q-chem gives an error
	for index in qm_atoms:
		#index+1 since q-chem starts with 1 as an index (compared to 0 in ase)
		f.write(str(index+1)+' ')
	f.write('\n$end\n\n')

def comments_section():
	'writes details of the comment section'
	g = open(create_input_dir+'/text-comments.txt','r')
	f.write(g.read()+'\n')
	g.close()

def ff_parameters():
	'writes details of the force field parameters'
	g = open(create_input_dir+'/text-ff.txt')
	f.write(g.read()+'\n')
	g.close()

def opt_section(fixed_atoms):
	'writes details of the $opt section'
	f.write('$opt\nfixed\n')
	for n in fixed_atoms:
		n = int(n)+1 #ase starts from zero but q-chem from 1
		f.write(str(n)+'\tXYZ\n')
	f.write('endfixed\n$end\n\n')

def molecules_section():
	'writes details of the #molecule section'
	f.write('$molecule\n0  '+str(multiplicity)+'\n')
	if os.path.getsize('tmp') == 0:
		print('Connectivity_NL results in an empty $mol section')
		exit()
	else:
		g = open('tmp', 'r')
		f.write(g.read())
		os.system('rm tmp')
		f.write('$end')

def qm_fixed_regions(traj, data, struc_dir):
	'qm region'
	qm_atoms    = data[traj+'.traj']['qm_region']
	fixed_atoms = []
	atoms       = io.read(struc_dir+'/'+traj+'.traj')
	for index, item in enumerate(atoms):
		if index in qm_atoms:
			continue
		else:
			fixed_atoms.append(str(index))
	return fixed_atoms, qm_atoms


def Al_Si_atoms():
	'finds number of Al/O/Si atoms in qm region'
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

def clean_input(input_file):
	'replaces the second line in input.xyz with a blank line'
	f = open(input_file,'r')
	lines = f.read().split("\n")
	f.close()

	g = open(input_file,'w')
	for index, line in enumerate(lines):
		if index == 1:
			g.write(' \n')
		elif line == '':
			continue
		else:
			g.write(line+"\n")
	g.close()

'Create dir for desired calculations'
if os.path.exists(calc_dir) is not True:
	os.system('mkdir '+calc_dir)

print('Calc dir', calc_dir)

for calc in calculations:
	
	print(calc)
	status = ''

	if job_type == 'sp':
		'check opt calc is done before generating sp folders'
		try:
			data_output[calc+'-opt-'+exchange+'-def2-svp-ref-'+data[calc+'.traj']['reference']]['energy']
		except:
			print(calc, 'opt calculation is incomplete')
			status = 'incomplete'
		
	if status == '':
		
		try:
			data[calc+'.traj']['reference']
		except:
			print('this is likely a duplicate that has been deleted')
			continue
		
			
		'avoid generating input for sp calculations where opt is incomplete'	
		if os.path.exists(calc_dir+'/'+calc+'-'+details+'-ref-'+data[calc+'.traj']['reference']) is not True:
			os.system('mkdir '+calc_dir+'/'+calc+'-'+details+'-ref-'+data[calc+'.traj']['reference'])
		os.chdir(calc_dir+'/'+calc+'-'+details+'-ref-'+data[calc+'.traj']['reference'])


		if basis == 'def2-sv(p)':
			os.system('cp '+struc_dir+'/'+calc+'.traj input.traj')
			atoms = io.read('input.traj')
			n_atoms = len(atoms)
			atoms.write('input.xyz')
			clean_input('input.xyz')
			
		elif basis == 'def2-tzvpd':
			os.system('cp '+calc_dir+'/'+calc+'-opt-'+exchange+'-def2-svp-ref-'+data[calc+'.traj']['reference']+'/full-atoms.xyz input.xyz')
			clean_input('input.xyz')

		try:
			os.system('cp '+data_dir+'/'+zeolite+'-NL.json .')
		except:
			print('original structure of zeolite-NL.json does not exist')
			exit()

		os.system('cp '+create_input_dir+'/connectivity_NL.py .')
		fixed_atoms, qm_atoms = qm_fixed_regions(calc, data, struc_dir)
		os.system('python connectivity_NL.py input.xyz '+str(qm_atoms)+' '+zeolite+' > tmp')

		with open("dir_data.json", "w") as write_file:
			json.dump(data[calc+'.traj'], write_file, indent=4)

		'''writing opt.in'''
		f = open('opt.in','w')
		rem_section()
		qm_atoms_section(qm_atoms)
		comments_section()
		ff_parameters()
		opt_section(fixed_atoms)
		molecules_section()
		f.close()

		'print initial qm structure'
		os.system('cp '+scripts_dir+'/qm_structure.py .')

		os.system('python qm_structure.py')	
		'add to dir_data number of QMMM atoms'	
		n_Al,n_Si, n_O = Al_Si_atoms()
		if n_Al == 1:
			H_qm = 12 #H added to qm region to justify MM region
		elif n_Al == 2:
			H_qm = H_qm_region(n_O, n_Si, n_Al)  #H added to qm region to justify MM region
		data[calc+'.traj']['QMMM length'] = len(data[calc+'.traj']['qm_region']) + H_qm #total atoms in qm region

		with open("dir_data.json", "w") as write_file:
			json.dump(data[calc+'.traj'], write_file, indent=4)


	exit()
