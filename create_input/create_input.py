#!/Users/hassanaljama/opt/anaconda3/bin/python

import os
from ase import io
import json

'''
Running calculations on selected strucutres. Provide Inputs below
'''
'Inputs'
calc         = [7,9,11]		#if continous list, input first and last. Otherwise, input individual entries
multiplicity = 2		#multiplicity of the structure
job_type     = 'opt' 		#either sp or opt
dir_name     = 'CHA-full-MR'	#name of the parent dir
zeolite      = 'CHA'		#zeolite name

'Inputs (rarely need a change)'
cwd              = os.getcwd()
exchange         = 'omegab97x-d'
struc_dir        = cwd+'/../structures_saved'
data_dir         = cwd+'/../data'
scripts_dir      = '/home/aljama/scripts/' 	#this is where qm-initial structure script is
calc_dir         = '/home/aljama/'+dir_name+'/calculations/'
create_input_dir = '/home/aljama/'+dir_name+'/create_input'

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
	elif job_type == 'sp':
		g = open(create_input_dir+'/text-rm-sp.txt','r')
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

'Create dir for desired calculations'

for calc in calculations:
	
	print(calc)
	status = ''

	if job_type == 'sp':
		'check opt calc is done before generating sp folders'
		try:
			data_output[calc+'-opt-omegab97x-d-def2-svp']['energy'] - 0
		except:
			print(calc, 'opt calculation is incomplete')
			status = 'incomplete'
		
	if status == '':
		'avoid generating input for sp calculations where opt is incomplete'	
		if os.path.exists(calc_dir+'/'+calc+'-'+details) is not True:
			os.system('mkdir '+calc_dir+'/'+calc+'-'+details)
		os.chdir(calc_dir+'/'+calc+'-'+details)

		with open("dir_data.json", "w") as write_file:
			json.dump(data[calc+'.traj'], write_file, indent=4)

		if basis == 'def2-sv(p)':
			os.system('cp '+struc_dir+'/'+calc+'.traj input.traj')
			atoms = io.read('input.traj')
			n_atoms = len(atoms)
			atoms.write('input.xyz')
		elif basis == 'def2-tzvpd':
			os.system('cp '+calc_dir+'/'+calc+'-opt-omegab97x-d-def2-svp/full-atoms.xyz input.xyz')

		try:
			os.system('cp '+data_dir+'/'+zeolite+'-NL.json .')
		except:
			print('original structure of zeolite-NL.json does not exist')
			exit()

		os.system('cp '+create_input_dir+'/connectivity_NL.py .')
		fixed_atoms, qm_atoms = qm_fixed_regions(calc, data, struc_dir)
		os.system('python connectivity_NL.py input.xyz '+str(qm_atoms)+' '+zeolite+' > tmp')

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

		'print structure of surrounding atoms'	
		os.system('cp '+scripts_dir+'/surroundings_structure.py .')
		os.system('python surroundings_structure.py')	
