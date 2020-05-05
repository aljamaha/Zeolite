import os

'Inputs'

run  = True
traj = False
cores = '8'
start = 1714
end   = start+4
calc_type = 'opt'
ref   = '93.traj'
cwd = '/home/aljama/BEA/H/calculations'

for i in range(start,end):
	i = str(i)
	if run == True:
		if calc_type == 'opt':
			os.chdir(cwd+'/'+i+'-opt-omegab97x-d-def2-svp-ref-'+ref)	
			os.system('run_qchem opt.in '+cores)
		elif calc_type == 'sp':
			try:
				os.chdir(cwd+'/'+i+'-sp-omegab97x-d-def2-tzvpd-ref'+ref)
				os.system('run_qchem opt.in '+cores)
			except:
				print(i+' does not exist')

	elif traj == True:
		if os.path.exists('traj'):
			os.system('cp '+cwd+'/'+i+'-opt-omegab97x-d-def2-svp-ref-'+ref+'/qm-final.xyz traj/'+i+'.xyz')	
		else:
			os.system('mkdir traj')
			os.system('cp '+cwd+'/'+i+'-opt-omegab97x-d-def2-svp-ref-'+ref+'/qm-final.xyz traj/'+i+'.xyz')	
		
	
