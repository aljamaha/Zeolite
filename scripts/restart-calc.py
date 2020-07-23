import os
os.chdir('/home/aljama/BEA/original/scripts')
from restart_opt import *

'Inputs'
exchange         = 'omegab97x-d'
basis            = 'def2-sv(p)'
multiplicity     = '1'
max_scf_cycles   = '150'
scf_algo         = 'diis_gdm'
calc_dir         = 'BEA/H'
THRESH		 = ''

create_input_dir = '/home/aljama/BEA/original/create_input/'
scripts_dir      = '/home/aljama/BEA/original/scripts/'
zeolite          = 'BEA'
run              = True

folders = ['1817-opt-omegab97x-d-def2-svp-ref-100.traj']

for f in folders:
	print(f)
	os.chdir('/home/aljama/'+calc_dir+'/calculations/'+f)
	os.system('ls | grep old')
	#os.system('tail opt.out')
	if run == True:
		try:
			os.system('mkdir old-1')
		except:
			pass
		os.system('cp *.* old-1')
		restart(exchange, basis, scf_algo, create_input_dir, scripts_dir, zeolite, multiplicity, max_scf_cycles, THRESH)
		#os.system('run_qchem opt.in 8')
