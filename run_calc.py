#!/Users/hassanaljama/opt/anaconda3/bin/python

import os
from ase import io

'''
Running calculations on selected strucutres. Provide Inputs below
'''

'Inputs'
calc = [111] 		#identify structures [under structures folder]
basis    = 'def2-sv(p)'	#basis set [def2-sv(p) or def2-tzvpd]
job_type = 'opt'
exchange = 'omegab97x-d'
fixed_atoms = '1400'
qm_atoms    = '14'

'General Inputs (do not change)'
cwd       = os.getcwd()
struc_dir = cwd+'/structures'
details  = job_type+'-'+exchange+'-'+basis.replace('(','').replace(')','') #naming dir (uniqueness)

def rm_section():
	'writes details of rm section'
	g = open(cwd+'/text-rm.txt','r')
	text_rm = g.read()
	g.close()
	f.write('$rm\n')
	f.write('jobtype\t'+job_type+'\n')
	f.write('exchange\t'+exchange+'\n')
	f.write('basis'+'\t'+basis+'\n')
	f.write('AIMD_FIXED_ATOMS\t'+fixed_atoms)
	f.write(text_rm)
	f.write('$end\n\n')

def qm_atoms_section(qm_atoms):
	'writes details of qm_atoms section'
	f.write('$qm_atoms\n1:'+qm_atoms+'\n$end\n\n')
	
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

def opt_section(qm_atoms, n_atom):
	'writes details of the $opt section'
	f.write('$opt\nfixed\n')
	for n in range(int(qm_atoms)+1,int(n_atoms)+1):
		f.write(str(n)+'\tXYZ\n')
	f.write('endfixed\n$end\n\n')
	
def molecules_section():
	'writes details of the #molecule section'
	f.write('$molecule\n0  1\n')
	g = open('tmp', 'r')
	f.write(g.read()+'\n')
	os.system('rm tmp')
	f.write('$end')

'run calculations'
for i in calc:
	i = str(i)
	os.chdir(cwd+'/calculations')
	if os.path.exists(i+'-'+details) is not True:
		os.system('mkdir '+i+'-'+details)
	os.chdir(i+'-'+details)
	os.system('cp '+struc_dir+'/'+i+'.traj input.traj')
	atoms = io.read('input.traj')
	n_atoms = len(atoms)
	atoms.write('input.xyz')
	os.system('cp '+cwd+'/qm-mm-connectivity.py .')
	os.system('python qm-mm-connectivity.py input.xyz '+qm_atoms+' > tmp')

	'''writing opt.in'''
	f = open('opt.in','w')

	rm_section()
	qm_atoms_section(qm_atoms)
	comments_section()
	ff_parameters()
	opt_section(int(qm_atoms), n_atoms)
	molecules_section()

	f.close()
