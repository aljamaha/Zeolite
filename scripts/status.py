import os, json, sys
from ase import io

'Gives a summary of the status of calculations'

'Inputs'
dir_name     = 'BEA/Pd1'

'Directroies'
calc_dir    = '/home/aljama/'+dir_name+'/calculations/'	#directory where caluculatiosn are saved
data_dir    = '/home/aljama/'+dir_name+'/data/'		#directory where data are saved

def folders_list(wd):
	'list of folders in a directory'
	os.chdir(wd)
	os.system("ls > tmp")
	folders = [line.rstrip('\n') for line in open('tmp')]
	folders.remove('tmp')
	os.system('rm tmp')

	return folders

def job_id_list(wd):
	'list of job ids under wd'
	try:
		os.chdir(wd)
		os.system("ls > tmp")
		folders = [line.rstrip('\n') for line in open('tmp')]
		folders.remove('tmp')
		os.system('rm tmp')
		job_id = []
		for name in folders:
			if len(name) == 5:
				job_id.append(name)
		return job_id

	except:
		pass

def running_jobs_list():
	'running jobs list'
	os.system('qstat -u aljama > tmp.txt')
	f = open('tmp.txt','r')
	jobs =  [line.rstrip('\n') for line in open('tmp.txt')]
	os.system('rm tmp.txt')
	running_job = []
	for job in jobs:
		for item in job.split(' '):
			if len(item) == 5:
				try:
					int(item)
					running_job.append(item)
				except:
					pass
	return running_job

def check_status(item):
	'check opt.out if calculation has failed'
	tmp = False
	lines =  [line.rstrip('\n') for line in open(calc_dir+'/'+item+'/opt.out')]
	for line in lines:
		if '**  MAXIMUM OPTIMIZATION CYCLES REACHED  **' in line:
			tmp = 'MAXIMUM'
			break
		elif 'Please submit' in line:
			tmp =  'FAILED'
			break
	return tmp

def comp(item):
	'check opt.out if calculation has completed'
	os.chdir(calc_dir+'/'+item)
	tmp = False
	try:
		os.system('tail opt.out > tmp')
		lines =  [line.rstrip('\n') for line in open(calc_dir+'/'+item+'/tmp')]
		os.system('rm tmp')
		for line in lines:
			if 'Thank you very much for using Q-Chem' in line:
				tmp = True
				break
	except:
		pass
	return tmp

def frozen_atoms():	
	'''
	checks if atoms are stuck (frozen)
	Outputs:
		True   - atoms are frozen
		False  - atoms are not frozen	 
	'''
	os.system('cat opt.out | grep NO > tmp')
	lines = [line.rstrip('\n') for line in open('tmp')]  #this is the best way
	k = 0
	for line in lines:
		if '***' in line:
			k+=1
	if k > 5:
		tmp = True
	else:
		tmp = False
	return tmp

folders = folders_list(calc_dir)	#folders in calculations/directory
running_jobs   = running_jobs_list()	#list of the running jobs
restart, failed, not_sure, frozen, Running, completed, terminate, terminate_id = [],[],[],[],[],[],[],[]

for item in folders:
	if 'def' in item:
		if os.path.exists(calc_dir+'/'+item+'/opt.out') == False:
			'Did not Start'
			print(item, '\t\t', 'Did not start') 
		elif comp(item) == True:
			'completed calculations'
			print(item, '\t\t', 'Completed')
			completed.append(item)
		else:
			try:	

				job_ids = job_id_list(calc_dir+'/'+item)
				mutual  = set(running_jobs) & set(job_ids)   

				if frozen_atoms() == True:
					'checks atoms are not stuck'
					if len(mutual) != 0:
						'Frozen and Running calculations'
						print(item, '\t\t', 'Must terminate')
						terminate.append(item)
						for i in mutual:
							terminate_id.append(int(i))
					else:
						'Frozen finished calculations'
						frozen.append(item)	
						print(item, '\t\t', 'Atoms are frozen. Adjust initial position')
				elif len(mutual) != 0:
					'Running calculations (not frozen)'
					print(item, '\t\t', 'Running')
					Running.append(item)	
				else:
					tmp_status = check_status(item)
					if tmp_status == 'MAXIMUM':
						'Reached Max optimization cycle'
						#if max_status() == 'pass':
						print(item, '\t\tMaximum opt reached. Restart')
						restart.append(item)
					elif tmp_status == 'FAILED':
						'Calc Failed'
						failed.append(item)
						print(item, '\t\tFailed')
					else:
						'Not sure!'
						print(item, 'Not sure!') 
						not_sure.append(item)
			except:
				print(item, '\t\tExcept failed')
				not_sure.append(item)

print('*******\nCalculations require restart:', len(restart), restart)
print('*******\nFrozen!:', len(frozen), frozen)
print('*******\nFailed calculations', len(failed), failed)
print('*******\nNot sure!:', len(not_sure), not_sure)
print('*******\nRunning:', len(Running), Running)
print('*******\nMust terminate:', len(terminate), terminate)
print('*******\nTerminate id:', terminate_id)
#print('Completed:', len(completed), completed)
