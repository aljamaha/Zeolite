#!/Users/hassanaljama/opt/anaconda3/bin/python

import os, json
from ase import io

'''
Create Q-Chem opt.in file
'''

class helpers():

	def __init__(self, exchange, basis, job_type, multiplicity, scf_algorithm, max_scf_cycles, wd, qm_atoms, traj, zeolite, cation, geom_opt_max_cycles, input_atoms, file_name, create_input_dir):

		'Main Inputs'
		self.exchange         = exchange
		self.basis	      = basis
		self.job_type         = job_type
		self.multiplicity     = multiplicity
		self.scf_algorithm    = scf_algorithm
		self.max_scf_cycles   = str(max_scf_cycles)
		self.qm_atoms	      = qm_atoms

		'optional'
		self.THRESH              = THRESH
		self.wd        	         = wd
		self.zeolite	         = zeolite
		self.cation              = cation
		self.traj	         = traj
		self.geom_opt_max_cycles = geom_opt_max_cycles
		self.file_name		 = file_name
		self.atoms 		 = input_atoms
		self.create_input_dir    = create_input_dir

		'defaults'
		if self.create_input_dir == '':
			self.create_input_dir = self.wd+'/'+self.zeolite+'/original/create_input'
		if self.file_name == '':
			self.file_name	      = 'opt.in'

		'qm atoms'
		if self.qm_atoms == '':
			'assumes qm atoms info in the json file'
			with open(self.wd+'/'+self.zeolite+'/'+self.cation+'/data/data.json', 'r') as read_file:
				data = json.load(read_file)
			self.qm_atoms    = data[self.traj+'.traj']['qm_region']

		'original structure'
		if self.atoms == '':
			self.atoms = io.read(self.wd+'/'+self.zeolite+'/original/structures_saved/'+self.traj+'.traj')
		else:
			self.atoms = io.read(self.atoms)

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
		g.write('geom_opt_coord\t0\n')
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
		if self.THRESH != '':
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
		g.write('endfixed\n$end\n')

	def molecules_section(self):
		'writes details of the #molecule section'
		g = open(self.file_name,'a')
		g.write('\n$molecule\n0  '+str(self.multiplicity)+'\n')
		if os.path.getsize('tmp') == 0:
			print('Connectivity_NL results in an empty $mol section')
			exit()
		else:
			f = open('tmp', 'r')
			f.write(g.read())
			os.system('rm tmp')
			f.write('$end')

def qchem_in(exchange, basis, multiplicity, job_type, wd='', THRESH='', qm_atoms ='', traj='', input_atoms='', zeolite='', geom_opt_max_cycles='', file_name='',create_input_dir):
	
	if wd == '':
		wd = os.getcwd()

	k = helpers(exchange, basis, job_type, multiplicity, scf_algorithm, max_scf_cycles, wd, qm_atoms, traj, zeolite, cation, geom_opt_max_cycles, input_atoms)

	'writing opt.in'
	f = open('opt.in','w')
	k.rem_section()
	k.qm_atoms_section()
	k.comments_section()
	k.ff_parameters_section()
	k.opt_section()
	f.close()

exchange       = 'omegab97x-d'	#'omegab97x-d' or 'B97-D3'
basis          = 'def2-sv(p)'	#'def2-sv(p) or def2-tzvpd'
multiplicity   = 2		#multiplicity of the structure
job_type       = 'sp' 		#either sp or opt
THRESH	       = '12'
scf_algorithm  = 'diis_gdm'
max_scf_cycles = '250'
wd	       = '/home/aljama/'
traj	       = '1'
zeolite        = 'BEA'		#zeolite name
cation         = 'Pd1'

qchem_in(exchange, basis, multiplicity, job_type, wd=wd, THRESH=THRESH, zeolite=zeolite, traj=traj)

