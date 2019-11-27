#!/Users/hassanaljama/opt/anaconda3/bin/python

import os

'''
Running calculations on selected strucutres
'''

'Inputs'
details  = 'B3LYP' 	#unique details about the calculations
basis    = '6-31+G**'	#basis set
job_type = 'opt'
exchange = 'omegab97x-d'
fixed_atoms = '1400'

'General Inputs'
cwd       = os.getcwd()
struc_dir = cwd+'/structures'

'identify structures to do calculations on'
calc = [111,112,113]

def rm_section():
	g = open(cwd+'/text-rm.txt','r')
	text_rm = g.read()
	g.close()
	f.write('$rm\n')
	f.write('jobtype\t'+job_type+'\n')
	f.write('exchange\t'+exchange+'\n')
	f.write('basis'+'\t'+basis+'\n$end\n')
	f.write('AIMD_FIXED_ATOMS\t'+fixed_atoms)
	f.write(text_rm)
	f.write('$end\n')

'run calculations'
for i in calc:
	i = str(i)
	os.chdir(cwd+'/calculations')
	if os.path.exists(i):
		continue
	else:
		os.system('mkdir '+i+'-'+details)
	os.chdir(i+'-'+details)
	os.system('cp '+struc_dir+'/'+i+'.traj input.traj')

	'''writing opt.in'''
	f = open('opt.in','w')
	rm_section()
	
	$qm_atoms
1:32
$end   



	f.close()
