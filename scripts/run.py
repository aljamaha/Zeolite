import os

'Inputs'

run  = False
traj = False
cores = ''
start = 
end   = 
calc_type = ''
ref   = '.traj'

cwd = os.getcwd()
os.system("ls > tmp")
folders = [line.rstrip('\n') for line in open('tmp')]  #this is the best way

for i in range(start,end):
	i = str(i)
	if run == True:
		if calc_type == 'opt':
			os.chdir(cwd+'/'+i+'-opt-omegab97x-d-def2-svp-ref-'+i+'.traj')	
			os.system('run_qchem opt.in '+cores)
		elif calc_type == 'sp':
			try:
				os.chdir(cwd+'/'+i+'-sp-omegab97x-d-def2-tzvpd')
				os.system('run_qchem opt.in '+cores)
			except:
				print(i+' does not exist')

	elif traj == True:
		if os.path.exists('traj'):
			os.system('cp '+cwd+'/'+i+'-opt-omegab97x-d-def2-svp-ref-'+ref+'/qm-final.xyz traj/'+i+'.xyz')	
		else:
			os.system('mkdir traj')
			os.system('cp '+cwd+'/'+i+'-opt-omegab97x-d-def2-svp-ref-'+ref+'/qm-final.xyz traj/'+i+'.xyz')	
		
	
