#!/Users/hassanaljama/opt/anaconda3/bin/python

import os, json
from ase import io

'''
Requirements:
		external_charge.py    : responsible for $external_charge section
		molecules_section.py  : reponseible for $molecules section
		sp calc. input.xyz - dir_data.json

		need to arrange fragments once calc. are done
'''

def restart(exchange, basis, scf_algo, scripts_dir, multiplicity, max_scf_cycles, THRESH=None):
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

	os.system('cp '+scripts_dir+'/convert-qchem-output-to-ase.py .')
	os.system('python convert-qchem-output-to-ase.py '+str(data["QMMM length"]))
	os.system('python external_charge.py > tmp')
	os.system('python molecules_section.py > tmp2')

	'''writing opt.in'''
	f = open('opt.in','w')
	rem_section(exchange, basis, multiplicity, max_scf_cycles, THRESH, scf_algo, f)
	molecules_section(multiplicity, f)
	external_charges(f)
	f.close()

def rem_section(exchange, basis, multiplicity, max_scf_cycles, THRESH, scf_algo, f):
	'writes details of rm section'
	f.write('$rem   \n')
	f.write('jobtype \teda\n')
	f.write('eda2\t2\n')
	f.write('exchange   '+exchange+'\n')
	f.write('basis   \t'+basis+'\n')
	f.write('max_scf_cycles  '+max_scf_cycles+'\n')         
	if THRESH != None:
		f.write('THRESH\t'+THRESH+'\n')
	f.write('scf_algorithm\t'+scf_algo+'\n')
	f.write('geom_opt_coords 0\n')
	f.write('geom_opt_max_cycles   	200\n')
	f.write('ecp	def2-ecp\n')
	f.write('symmetry    	off\n')
	f.write('sym_ignore   	true\n')
	f.write('print_input   	true\n')
	f.write('mem_total   	8000\n')
	f.write('mem_static   	450\n')
	f.write('geom_opt_dmax   80\n')
	f.write('pop_mulliken false\n')
	f.write('scf_print_frgm true\n')
	f.write('$end\n\n')

def molecules_section(multiplicity, f):
	'writes details of the #molecule section'
	f.write('$molecule\n0  '+multiplicity+'\n')
	f.write('--\n')
	g = open('tmp2', 'r')
	f.write(g.read())
	os.system('rm tmp2')
	f.write('$end\n\n')

def external_charges(f):
	'writes details of the #molecule section'
	f.write('$external_charges\n')
	g = open('tmp', 'r')
	f.write(g.read())
	os.system('rm tmp')
	f.write('$end')

exchange 	 = 'omegab97x-d'
basis    	 = 'def2-tzvpd'
scf_algo	 = 'diis_gdm'
scripts_dir      = '/home/aljama/BEA/original/scripts/'
multiplicity     = '1'
max_scf_cycles   = '150'
THRESH		 = '12'
restart(exchange, basis, scf_algo,  scripts_dir, multiplicity, max_scf_cycles, THRESH=THRESH)
