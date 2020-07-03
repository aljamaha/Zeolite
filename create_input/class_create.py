#!/Users/hassanaljama/opt/anaconda3/bin/python

import os
from ase import io
import json

'''
Running calculations on selected strucutres. Provide Inputs below
'''

class create_input():

	#def __init__(self, dir_name, multiplicity, exchange, job_type, zeolite, THRESH, THRESH_value, scf_algorithm, calc,wd):
	def __init__(self, dir_name, exchange, job_type):
		'Main Inputs'
		self.dir_name         = dir_name
		#self.multiplicity     = multiplicity
		self.exchange         = exchange
		self.job_type         = job_type
		#self.zeolite          = zeolite
		#self.THRESH           = THRESH
		#self.THRESH_value     = THRESH_value
		#self.scf_algorithm    = scf_algorithm
		#self.calc	      = calc
		#self.wd		      = wd
		#self.struc_dir        = wd+'/'+dir_name+'/structures_saved'
		self.data_dir	      = wd+'/'+dir_name+'/data'
		#self.scripts_dir      = wd+'/scripts'
		#self.calc_dir         = wd+'/'+dir_name+'/calculations'
		#self.create_input_dir = wd+'/'+zeolite+'/original/create_input'
		
	def basis_set(self):
		'basis set depending on job_type'
		if self.job_type == 'sp':
			basis    = 'def2-tzvpd'
		elif self.job_type == 'opt':
			basis    = 'def2-sv(p)'
		details  = self.job_type+'-'+self.exchange+'-'+basis.replace('(','').replace(')','') #naming dir (uniqueness)
		return basis

	def calc_list(self, calc):
		'''
		Indicies of calculations to be conducted'
		Input: calc - calculations indicies (either start and end point, or a complete list)
		Output: complete list
		'''
		calculations = []
		if len(calc) == 2:
			'continous list from start to end'
			for j in range(calc[0],calc[1]):
				calculations.append(str(j))
		else:
			'individual list'
			for j in calc:
				calculations.append(str(j))
		return calculations

	def load_data(self):
		'Load json data (data.json and data_output.json)'
		with open(self.data_dir+"/data.json", "r") as read_file:
		    data = json.load(read_file)
		try:
			with open(self.data_dir+"/data_output.json", "r") as read_file:
				data_output = json.load(read_file)
		except:
			pass

		return data, data_output

''' Testing ... '''
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
k = create_input(dir_name, exchange, job_type)
print(k.basis_set())
print(k.calc_list(calc))
k.load_data()

'''
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
	f.write('scf_algorithm\t'+scf_algorithm+'\n')
	f.write(text_rm)
	f.write('model_system_mult '+str(multiplicity)+'\n')
	if THRESH == True:
		f.write('THRESH\t'+str(THRESH_value)+'\n')
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
	#''identifies the number of H atoms in qm region added to justify dangling bonds
	Inputs:
		n_O:  number of O  in qm region
		n_Si: number of Si in qm region
		n_Al: number of Al in qm region
	Outputs:
		number of additional H in qm region to prevent dangling Si bonds
	#''

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

			if exchange == 'omegab97x-d':	
				try:
					'copy GGA results to hGGA if available'	
					wd = os.getcwd()
					os.chdir(calc_dir+'/'+calc+'-opt-B97-D3-def2-svp-ref-'+data[calc+'.traj']['reference'])
					os.system('cp '+scripts_dir+'/qchem-to-ase-all-atoms.py .')
					os.system('python qchem-to-ase-all-atoms.py')
					os.system('rm -r traj-full-atoms/')
					os.chdir(wd)
					os.system('cp '+calc_dir+'/'+calc+'-opt-B97-D3-def2-svp-ref-'+data[calc+'.traj']['reference']+'/full-atoms.xyz input.xyz')
					atoms = io.read('input.xyz')
					print('optimized GGA calc exists')
				except:
					print('No optimized GGA calc')
					atoms = io.read('input.traj')
					atoms.write('input.xyz')
			else:	
				atoms = io.read('input.traj')
				atoms.write('input.xyz')

			n_atoms = len(atoms)
			clean_input('input.xyz')
			
		elif basis == 'def2-tzvpd':
			wd = os.getcwd()
			os.chdir(calc_dir+'/'+calc+'-opt-'+exchange+'-def2-svp-ref-'+data[calc+'.traj']['reference'])
			os.system('cp '+scripts_dir+'/qchem-to-ase-all-atoms.py .')
			os.system('python qchem-to-ase-all-atoms.py')
			os.system('rm -r traj-full-atoms/')
			os.chdir(wd)
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

		#writing opt.in
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

'''
