#!/Users/hassanaljama/opt/anaconda3/bin/python

import os, json
from ase import io

def restart(exchange, basis, scf_algo, create_input_dir, scripts_dir, Zeolite, multiplicity, max_scf_cycles, THRESH=''):
	'''
	Generates a new opt.in based on terminated opt.out calculation
	* assumes the calculation to be restarted is in the same folder
	Inputs:
		exchange
		basis set
		scf_algo
		create input dir
		scripts dir
		zeolite
		multiplicity
		max_scf_cycles
	'''

	'Load information on qm region'
	with open("dir_data.json", "r") as read_file:
		data = json.load(read_file)

	'Run the code'
	atoms = io.read('input.xyz')
	n_atoms = len(atoms)
	os.system('cp '+scripts_dir+'/qchem-to-ase-all-atoms.py .')
	os.system('python qchem-to-ase-all-atoms.py opt.out input.xyz '+str(n_atoms))
	os.system('cp '+create_input_dir+'/connectivity_NL.py .')
	fixed_atoms, qm_atoms = qm_fixed_regions(data, atoms)
	os.system('python connectivity_NL.py full-atoms.xyz '+str(qm_atoms)+' '+Zeolite+' > tmp')

	'''writing opt.in'''
	f = open('opt.in','w')
	rem_section(create_input_dir, exchange, basis, fixed_atoms, multiplicity, max_scf_cycles, THRESH, scf_algo, f)
	qm_atoms_section(qm_atoms, f)
	comments_section(create_input_dir, f)
	ff_parameters(create_input_dir, f)
	opt_section(fixed_atoms, f)
	molecules_section(multiplicity, f)
	f.close()

def rem_section(create_input_dir, exchange, basis, fixed_atoms, multiplicity, max_scf_cycles, THRESH, scf_algo, f):
	'writes details of rm section'
	g = open(create_input_dir+'/text-rm.txt','r')
	text_rm = g.read()
	g.close()
	f.write('$rem   \n')
	f.write('jobtype \topt\n')
	f.write('exchange   '+exchange+'\n')
	f.write('basis   \t'+basis+'\n')
	f.write('AIMD_FIXED_ATOMS \t'+str(len(fixed_atoms))+'\n')
	f.write('model_system_mult '+multiplicity+'\n')
	f.write('max_scf_cycles  '+max_scf_cycles+'\n')         
	if THRESH == True:
		f.write('THRESH\t14\n')
	f.write('scf_algorithm\t'+scf_algo+'\n')
	f.write(text_rm)
	f.write('$end\n\n')

def qm_atoms_section(qm_atoms, f):
	'writes details of qm_atoms section'
	f.write('$qm_atoms\n')
	qm_atoms.sort() #if unsorted q-chem gives an error
	for index in qm_atoms:
		#index+1 since q-chem starts with 1 as an index (compared to 0 in ase)
		f.write(str(index+1)+' ')
	f.write('\n$end\n\n')

def comments_section(create_input_dir, f):
	'writes details of the comment section'
	g = open(create_input_dir+'/text-comments.txt','r')
	f.write(g.read()+'\n')
	g.close()

def ff_parameters(create_input_dir, f):
	'writes details of the force field parameters'
	g = open(create_input_dir+'/text-ff.txt')
	f.write(g.read()+'\n')
	g.close()

def opt_section(fixed_atoms, f):
	'writes details of the $opt section'
	f.write('$opt\nfixed\n')
	for n in fixed_atoms:
		n = int(n)+1 #ase starts from zero but q-chem from 1
		f.write(str(n)+'\tXYZ\n')
	f.write('endfixed\n$end\n\n')

def molecules_section(multiplicity, f):
	'writes details of the #molecule section'
	f.write('$molecule\n0  '+multiplicity+'\n')
	g = open('tmp', 'r')
	f.write(g.read())
	os.system('rm tmp')
	f.write('$end')

def qm_fixed_regions(data, atoms):
	'qm region'
	qm_atoms    = data['qm_region']
	fixed_atoms = []
	for index, item in enumerate(atoms):
		if index in qm_atoms:
			continue
		else:
			fixed_atoms.append(str(index))
	return fixed_atoms, qm_atoms
