import os, json, sys
from ase import io

'Gives a summary of the status of calculations'

'Inputs'
dir_name     = 'BEA/H'

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

def max_status():
	'''
	checks if calculation that run max steps is stuck (atoms are freezing)
	Outputs:
		failed - atoms are frozen
		pass   - need a restart
	'''
	os.system('cat opt.out | grep NO > tmp')
	lines = [line.rstrip('\n') for line in open('tmp')]  #this is the best way
	k = 0
	for line in lines:
		if '***' in line:
			k+=1
	if k > 5:
		tmp = 'fail'
	else:
		tmp = 'pass'
	return tmp

'Load data'
try:
	with open(data_dir+"/data_output.json", "r") as read_file:
		data = json.load(read_file)		#output data
except:
	print('Data not available')

folders = folders_list(calc_dir)	#folders in calculations/directory
running_jobs   = running_jobs_list()	#list of the running jobs
restart = []
for item in data:
	if 'def' in item:
		if data[item]['status'] == 'complete':
			'completed calculations'
			print(item,'\t\t', 'Completed')
		else:
			job_ids = job_id_list(calc_dir+'/'+item)
			try:
				mutual  = set(running_jobs) & set(job_ids)   
				if len(mutual) != 0:
					'Running calculations'
					print(item, '\t\t', 'Running')	
				else:
					if os.path.exists(calc_dir+'/'+item+'/opt.out') == False:
						'Did not start Calculations'
						print(item,'\t\t', 'Did not start')
					else:
						if check_status(item) == 'MAXIMUM':
							'Reached Max optimization cycle'
							if max_status() == 'pass':
								print(item, '\t\t', 'Maximum opt cycle is reached. Need a restart')
								restart.append(item)
							else:
								print(item, '\t\t', 'Atoms are frozen. Adjust initial position')
						elif check_status(item) == 'FAILED':
							'Calc Failed'
							print(item, '\t\t', 'Failed')
						else:
							'Not sure!'
							print(item, '????') 
			except:
				pass

print('Calculations require restart:', restart)
