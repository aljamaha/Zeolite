#!/Users/hassanaljama/opt/anaconda3/bin/python

import os, json
from ase import io

'''
Create Q-Chem opt.in file
'''

class helpers():

	def __init__(self, dir_name, multiplicity, exchange, job_type, zeolite, THRESH, scf_algorithm, calc,wd, max_scf_cycles, traj):
		'Main Inputs'
		self.dir_name         = dir_name
		self.multiplicity     = multiplicity
		self.exchange         = exchange
		self.job_type         = job_type
		self.zeolite          = zeolite
		self.THRESH           = THRESH
		self.scf_algorithm    = scf_algorithm
		self.calc	      = calc
		self.wd		      = wd
		self.struc_dir        = wd+'/'+dir_name+'/structures_saved'
		self.data_dir	      = wd+'/'+dir_name+'/data'
		self.scripts_dir      = wd+'/scripts'
		self.calc_dir         = wd+'/'+dir_name+'/calculations'
		self.create_input_dir = wd+'/'+zeolite+'/original/create_input'
		self.file_name	     = 'opt.in'
		self.max_scf_cycles = str(max_scf_cycles)
		self.traj	 = traj

		with open(self.data_dir+"/data.json", "r") as read_file:
			self.data = json.load(read_file)
		self.qm_atoms    = self.data[self.traj+'.traj']['qm_region']
		self.fixed       = self.fixed_atoms()
		print(len(self.fixed))

	'''
	def fixed_atoms(self):
		fixed = []
		atoms       = io.read(self.struc_dir+'/'+traj+'.traj')
		for index, item in enumerate(atoms):
			if index not in self.qm_atoms:
				fixed.append(str(index))
		return fixed
	

	#def data_output(sef):
	#	try:
	#		with open(self.data_dir+"/data_output.json", "r") as read_file:
	#			data_output = json.load(read_file)
	#			return data_output
	#	except:
	#		pass
	'''

	def rem_section(self):
		'writes details of rm section'
		g = open(self.file_name,'a')
		g.write('$rem   \n')
		g.write('geom_opt_coord\t0\n')
		g.write('max_scf_cycles\t'+self.max_scf_cycles+'\n')
		g.write('geom_opt_max_cycles\t150\n')
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
		g.write('jobtype \t'+self.job_type+'\n')
		g.write('exchange   '+self.exchange+'\n')
		g.write('scf_algorithm\t'+self.scf_algorithm+'\n')
		g.write('model_system_mult '+str(self.multiplicity)+'\n')
		if self.THRESH != '':
			g.write('THRESH\t'+str(self.THRESH)+'\n')
		#self.fixed = self.fixed_atoms()
		g.write('AIMD_FIXED_ATOMS \t'+str(len(self.fixed))+'\n')
		#g.write('basis   \t'+basis+'\n')
		g.write('$end\n\n')
		g.close()

	def comments(self):
		'writes comments related to QM/MM'
		g = open(self.file_name,'a')
		g.write('\n')
		lines = [line for line in open(self.create_input_dir+'/text-comments.txt')]
		for line in lines:
			g.write(line)
		g.write('\n')
		g.close()

	def ff_parameters(self):
		'writes details of the force field parameters'
		g = open(self.file_name,'a')
		g.write('\n')
		lines = [line for line in open(self.create_input_dir+'/text-ff.txt')]
		for line in lines:
			g.write(line)
		g.write('\n')
		g.close()

	def fixed_atoms(self):
		'fixed region'
		#self.fixed = []
		fixed = []
		atoms = io.read(self.struc_dir+'/'+self.traj+'.traj')
		for index, item in enumerate(atoms):
			if index in self.qm_atoms:
				continue
			else:
				#self.fixed.append(str(index))
				fixed.append(str(index))
		
		return fixed

def qchem_in():
	wd	      = '/home/aljama/'
	dir_name      = 'BEA/Pd1'	#name of the parent dir
	multiplicity  = 2		#multiplicity of the structure
	exchange      = 'omegab97x-d'	#'omegab97x-d' or 'B97-D3'
	job_type      = 'sp' 		#either sp or opt
	zeolite       = 'BEA'		#zeolite name
	THRESH	      = False
	THRESH_value  = ''
	scf_algorithm = 'diis_gdm'
	calc	      = [1309]
	max_scf_cycles = '250'
	traj	       = '1'

	k = helpers( dir_name, multiplicity, exchange, job_type, zeolite, THRESH, scf_algorithm, calc,wd, max_scf_cycles, traj)
	'writing opt.in'
	f = open('opt.in','w')
	k.rem_section()
	k.comments()
	k.ff_parameters()
	f.close()

qchem_in()

