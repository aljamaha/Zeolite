#!/Users/hassanaljama/opt/anaconda3/bin/python

import os, json
from ase import io
from functions import *
from connectivity import connect

'''
Create Q-Chem opt.in file
'''

class input_text():

	def __init__(self, exchange, basis, job_type, multiplicity, scf_algorithm, create_input_dir, max_scf_cycles, wd, THRESH, qm_atoms, traj, zeolite, cation, geom_opt_max_cycles, initial_structure, file_name):

		'Main Inputs'
		self.exchange            = exchange		#exchange function
		self.basis	         = basis			#basis set
		self.job_type            = job_type		#job type (opt, sp, freq, etc.)
		self.multiplicity        = multiplicity		#multiplicity
		self.scf_algorithm       = scf_algorithm		#diis or diis_gdm
		self.qm_atoms	         = qm_atoms		#atoms in the qm region
		self.create_input_dir    = create_input_dir	#directory where .txt files are present
		self.initial_structure 	 = initial_structure	#name of .xyz file containing reactant (if using database, use input.xyz)

		'Required if using database'
		self.zeolite	         = zeolite		#type of zeolite material
		self.cation              = cation		#name of the cation
		self.wd        	         = wd			#working directory
		self.traj	         = traj			#item number of traj file in database

		'optional'
		self.THRESH              = THRESH		#THRESH value (12-15)
		self.geom_opt_max_cycles = geom_opt_max_cycles	#150-300
		self.max_scf_cycles      = str(max_scf_cycles)	#max scf cycles (150-300)
		self.file_name		 = file_name		#name of the file (default is opt.in)

		'qm atoms'
		if self.qm_atoms == None:
			try:
				'assumes qm atoms info in the json file based on default arrangement of the code'
				with open(self.wd+'/'+self.zeolite+'/'+self.cation+'/data/data.json', 'r') as read_file:
					data = json.load(read_file)
				self.qm_atoms    = data[self.traj+'.traj']['qm_region']
			except:
				print('qm atoms could not be found')

		'original structure'
		if self.initial_structure == 'input.xyz' and cation != None:
			self.atoms = io.read(self.wd+'/'+self.zeolite+'/'+cation+'/structures_saved/'+self.traj+'.traj')
		else:
			self.atoms = io.read(self.initial_structure)

	def fixed_atoms(self):
		'fixed atoms'
		fixed = []
		for index, item in enumerate(self.atoms):
			if index not in self.qm_atoms:
				fixed.append(str(index))
		return fixed

	def rem_section(self):
		'writes details of rm section'
		g = open(self.file_name,'a')
		g.write('$rem   \n')
		g.write('jobtype \t'+self.job_type+'\n')
		g.write('exchange   '+self.exchange+'\n')
		g.write('basis   \t'+self.basis+'\n')
		g.write('model_system_mult '+str(self.multiplicity)+'\n')
		g.write('max_scf_cycles\t'+self.max_scf_cycles+'\n')
		g.write('scf_algorithm\t'+self.scf_algorithm+'\n')
		g.write('AIMD_FIXED_ATOMS \t'+str(len(self.fixed_atoms()))+'\n')
		g.write('geom_opt_coords\t0\n')
		g.write('geom_opt_max_cycles\t'+str(self.geom_opt_max_cycles)+'\n')
		g.write('ecp\tdef2-ecp\n')
		g.write('QM_MM_INTERFACE Zeolite\n')
		g.write('force_field   	charmm27\n')
		g.write('user_connect   	true\n')
		g.write('symmetry    	off\n')
		g.write('sym_ignore   	true\n')
		g.write('print_input   	true\n')
		g.write('qmmm_print   	false\n')
		g.write('qm_mm   	true\n')
		g.write('qmmm_full_hessian   false\n')
		g.write('mem_total   	8000\n')
		g.write('mem_static   	450\n')
		g.write('geom_opt_dmax   80\n')
		g.write('pop_mulliken false\n')
		if self.THRESH != None:
			g.write('THRESH\t'+str(self.THRESH)+'\n')
		g.write('$end\n')
		g.close()

	def qm_atoms_section(self):
		'writes details of qm_atoms section'
		g = open(self.file_name,'a')
		g.write('\n$qm_atoms\n')
		self.qm_atoms.sort() #if unsorted q-chem gives an error
		for index in self.qm_atoms:
			#index+1 since q-chem starts with 1 as an index (compared to 0 in ase)
			g.write(str(index+1)+' ')
		g.write('\n$end\n')
		g.close()

	def comments_section(self):
		'writes comments related to QM/MM'
		g = open(self.file_name,'a')
		g.write('\n')
		try:
			lines = [line for line in open(self.create_input_dir+'/text-comments.txt')]
		except:
			print('Could not find text-comments.txt file in create input dir')
			exit()
		for line in lines:
			g.write(line)
		g.write('\n')
		g.close()

	def ff_parameters_section(self):
		'writes details of the force field parameters'
		g = open(self.file_name,'a')
		g.write('\n')
		lines = [line for line in open(self.create_input_dir+'/text-ff.txt')]
		for line in lines:
			g.write(line)
		g.write('\n')
		g.close()

	def opt_section(self):
		'writes details of the $opt section'
		g = open(self.file_name,'a')
		g.write('\n$opt\nfixed\n')
		for n in self.fixed_atoms():
			n = int(n)+1 #ase starts from zero but q-chem from 1
			g.write(str(n)+'\tXYZ\n')
		g.write('endfixed\n$end\n\n')

	def molecules_section(self):
		g = open(self.file_name,'a')
		'writes details of the #molecule section'
		g.write('$molecule\n0  '+str(self.multiplicity)+'\n')
		try:
			'copying NL dir'
			os.system('cp '+self.wd+'/'+self.zeolite+'/'+self.cation+'/data/*-NL.json .')
		except:
			pass

		if self.zeolite != None:
			connect(self.initial_structure, self.qm_atoms, zeolite=self.zeolite)
		else:
			connect(self.initial_structure, self.qm_atoms)
		
		if os.path.getsize('tmp') == 0:
			print('Connectivity_NL results in an empty $mol section')
			exit()
		else:
			r = open('tmp', 'r')
			g.write(r.read())
			os.system('rm tmp')
			g.write('$end')
			r.close()

def qchem_input(exchange, basis, job_type, multiplicity, scf_algorithm, create_input_dir, initial_structure='input.xyz', cation=None, max_scf_cycles='150',  wd=None, THRESH=None, qm_atoms=None, traj=None, zeolite=None, geom_opt_max_cycles=150, file_name='opt.in'):

	'''
	creats input file for qchem calculations
	
	Required inputs:
	exchange		#exchange function
	basis			#basis set
	job_type		#job type (opt, sp, freq, etc.)
	multiplicity		#multiplicity
	scf_algorithm		#diis or diis_gdm
	qm_atoms		#list of atoms in the qm region
	create_input_dir	#directory where .txt files are present (rem.txt etc.)

	Optional inputs:
	THRESH			#THRESH value (8-15)
	str(max_scf_cycles)	#max scf cycles (150-300)
	wd			#working directory
	zeolite			#type of zeolite material
	cation			#name of the cation
	traj			#item number of in the json database
	geom_opt_max_cycles	#150-300
	file_name		#name of the file (default is opt.in)
	initial_structure	#name of .xyz file containing reactant
	'''
	
	k = input_text(exchange, basis, job_type, multiplicity, scf_algorithm,  create_input_dir, max_scf_cycles, wd, THRESH, qm_atoms, traj, zeolite, cation, geom_opt_max_cycles, initial_structure, file_name)

	'writing opt.in'
	f = open('opt.in','w')
	k.rem_section()
	k.qm_atoms_section()
	k.comments_section()
	k.ff_parameters_section()
	k.opt_section()
	k.molecules_section()
	f.close()

'''
'Test it here .. '
exchange       = 'B97-D3'	#'omegab97x-d' or 'B97-D3'
basis          = 'def2-sv(p)'	#'def2-sv(p) or def2-tzvpd'
multiplicity   = 2		#multiplicity of the structure
job_type       = 'opt' 		#either sp or opt
THRESH	       = ''
scf_algorithm  = 'diis'
max_scf_cycles = '250'
wd	       = '/home/aljama/'
traj	       = '2055'
zeolite        = 'BEA'		#zeolite name
#cation         = 'Pd1'

create_input_dir =  '/home/aljama/BEA/original/create_input/'
with open('dir_data.json','r') as read:
	data = json.load(read)
	data = data['qm_region']

qchem_input(exchange, basis, job_type, multiplicity, scf_algorithm, create_input_dir, initial_structure='a.xyz', cation='Pd1',wd=wd, traj='1100',zeolite='BEA')
'''
