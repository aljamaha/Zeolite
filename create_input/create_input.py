#!/Users/hassanaljama/opt/anaconda3/bin/python

import os
from ase import io
import json

'''
Running calculations on selected strucutres. Provide Inputs below
'''

'Inputs'
calc     = ['1']	#identify structures [under structures folder]
basis    = 'def2-sv(p)'	#basis set [def2-sv(p) or def2-tzvpd]
job_type = 'opt'
exchange = 'omegab97x-d'

'General Inputs (do not change)'
cwd       = os.getcwd()
struc_dir = cwd+'/../structures'
details  = job_type+'-'+exchange+'-'+basis.replace('(','').replace(')','') #naming dir (uniqueness)
data_dir = cwd+'/../data'
with open(data_dir+"/data.json", "r") as read_file:
    data = json.load(read_file)
		
def rm_section():
	'writes details of rm section'
	g = open(cwd+'/text-rm.txt','r')
	text_rm = g.read()
	g.close()
	f.write('$rm\n')
	f.write('jobtype\t'+job_type+'\n')
	f.write('exchange\t'+exchange+'\n')
	f.write('basis'+'\t'+basis+'\n')
	f.write('AIMD_FIXED_ATOMS\t'+str(len(fixed_atoms))+'\n')
	f.write(text_rm)
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
	g = open(cwd+'/text-comments.txt','r')
	f.write(g.read()+'\n')
	g.close()

def ff_parameters():
	'writes details of the force field parameters'
	g = open(cwd+'/text-ff.txt')
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
	f.write('$molecule\n0  1\n')
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

'run calculations'

for i in calc:
	i = str(i)
	os.chdir(cwd+'/../calculations')
	if os.path.exists(i+'-'+details) is not True:
		os.system('mkdir '+i+'-'+details)
	os.chdir(i+'-'+details)
	os.system('cp '+struc_dir+'/'+i+'.traj input.traj')
	atoms = io.read('input.traj')
	n_atoms = len(atoms)
	atoms.write('input.xyz')
	os.system('cp '+cwd+'/../general/connectivity.py .')
	fixed_atoms, qm_atoms = qm_fixed_regions(i, data, struc_dir)
	os.system('python connectivity.py input.xyz '+str(qm_atoms)+' > tmp')

	'''writing opt.in'''
	f = open('opt.in','w')
	rm_section()
	qm_atoms_section(qm_atoms)
	comments_section()
	ff_parameters()
	opt_section(fixed_atoms)
	molecules_section()
	f.close()
