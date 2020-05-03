#!/Users/hassanaljama/opt/anaconda3/bin/python

import os
from ase import io
import json

'''
create input files with only qm region (to compare nuclear repulsion energies)
'''

'Inputs'
calc         = [1,2000]	#if continous list, input first and last. Otherwise, input individual entries
multiplicity = 1	#multiplicity of the structure
job_type     = 'opt' 	#either sp or opt
dir_name     = 'BEA/BEA-full-Si'	#name of the parent dir
zeolite      = 'BEA'	#zeolite name

'Inputs (rarely need a change)'
cwd              = os.getcwd()
exchange         = 'omegab97x-d'
struc_dir        = cwd+'/../structures_saved'
data_dir         = cwd+'/../data'
scripts_dir      = '/home/aljama/scripts/' 	#this is where qm-initial structure script is
calc_dir         = '/home/aljama/'+dir_name+'/calculations/'
create_input_dir = '/home/aljama/'+dir_name+'/create_input/'

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

def rem_section():
	'writes details of rm section'
	if job_type == 'opt':
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

def molecules_section(charge):
	'writes details of the #molecule section'
	f.write('$molecule\n'+str(charge)+'  '+str(multiplicity)+'\n')
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
	fixed_atoms, qm_atoms = [],[]
	for atom in atoms:
		qm_atoms.append(atom.index)	
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

for calc in calculations:
	
		print(calc)
		
		'check calc is available'
		try:
			data[calc+'.traj']['reference']
		except:
			print('this is likely a duplicate that has been deleted')
			continue
			
		'create folder'
		if os.path.exists(calc_dir+'/'+calc+'-'+details+'-ref-'+data[calc+'.traj']['reference']) is not True:
			os.system('mkdir '+calc_dir+'/'+calc+'-'+details+'-ref-'+data[calc+'.traj']['reference'])
		os.chdir(calc_dir+'/'+calc+'-'+details+'-ref-'+data[calc+'.traj']['reference'])

		'copy main input traj'
		os.system('cp '+struc_dir+'/'+calc+'.traj input.traj')
		atoms = io.read('input.traj')
		atoms.write('input.xyz')

		try:
			os.system('cp '+data_dir+'/'+zeolite+'-NL.json .')
		except:
			print('original structure of zeolite-NL.json does not exist')

		'print initial qm structure'
		with open("dir_data.json", "w") as write_file:
			json.dump(data[calc+'.traj'], write_file, indent=4)
		os.system('cp '+scripts_dir+'/qm_structure.py .')
		os.system('python qm_structure.py')	
		atoms = io.read('qm-initial.traj')
		atoms.write('input.xyz')

		'connectivity'
		os.system('cp '+create_input_dir+'/connectivity_NL.py .')
		fixed_atoms, qm_atoms = qm_fixed_regions(calc, data, struc_dir)
		n_atoms = len(atoms)
		os.system('python connectivity_NL.py input.xyz '+str(qm_atoms)+' '+zeolite+' > tmp')

		'charge'
		n_Al,n_Si, n_O = Al_Si_atoms()
		if n_Al == 1:
			charge = '1'	
		elif n_Al == 2:
			charge = '2'

		'''writing opt.in'''
		f = open('opt.in','w')
		rem_section()
		qm_atoms_section(qm_atoms)
		ff_parameters()
		opt_section(fixed_atoms)
		molecules_section(charge)
		f.close()

		with open("dir_data.json", "w") as write_file:
			json.dump(data[calc+'.traj'], write_file, indent=4)
