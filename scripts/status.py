import os, json, sys
from ase import io
from copy import deepcopy

'Gives a summary of the status of calculations'

'Inputs'
try:
	zeolite  = sys.argv[1]
	ads      = sys.argv[2]
except:
	print('python status.py [Zeolite] [ads]')
	exit()

'Directroies'
calc_dir    = '/home/aljama/'+zeolite+'/'+ads+'/calculations/'	#directory where caluculatiosn are saved
data_dir    = '/home/aljama/'+zeolite+'/'+ads+'/data/'		#directory where data are saved
cwd = os.getcwd()

data = {}

def MR4(cwd):
	'opt.out for 4MR have an added line of 4MR at the end'
	tmp = False
	os.chdir(cwd)
	os.system('tail opt.out | grep 4MR > tmp_4MR')
	if os.path.getsize('tmp_4MR') != 0:
		tmp = True
	return tmp

def diverged(cwd):
	'determines if the calculation diverged'

	status = True #assumes it did not diverge

	os.chdir(cwd)
	os.system('cat opt.out | grep "Energy change" > tmp')
	lines = [line.rstrip('\n') for line in open('tmp')]
	E, failures = [],0

	for index, line in enumerate(lines):

		try:
			'extract energy change per line'
			tmp = float(line[25:-20])
			E.append(tmp)
		except:
			tmp = False

		if tmp != False and index not in [0,1]:
			'if there is an energy in the line (and not ***)'
			diff = E[-1] - E[-2]
			if diff < 0 and abs(diff) > 1e-5:
				failures += 1

	if failures > 20:
		status = False

	return status

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
	os.system('tail '+calc_dir+'/'+item+'/opt.out > tmp')
	lines =  [line.rstrip('\n') for line in open(calc_dir+'/'+item+'/tmp')]
	os.system('rm tmp')
	
	for line in lines:
		if '**  MAXIMUM OPTIMIZATION CYCLES REACHED  **' in line:
			tmp = 'MAXIMUM'
			break
		elif 'SCF failed to converge' in line:
			tmp = 'SCF failed to converge' 
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
restart, failed, not_sure, frozen, Running, completed, terminate, terminate_id = [],[],[],[],[],[],[],''
scf_converge, not_started, diverge, four_MR = [],[],[],[]

for item in folders:
	if 'def' in item:
		if os.path.exists(calc_dir+'/'+item+'/opt.out') == False:
			'Did not Start'
			#print(item, '\t\t', 'Did not start') 
			not_started.append(item)
			data[item] = 'Did not start'	
		elif MR4(calc_dir+'/'+item) == True:
			'initial condition of cation placed in 4 MR'
			four_MR.append(item)
			data[item] = '4MR'
		elif comp(item) == True:
			'completed calculations'
			#print(item, '\t\t', 'Completed')
			completed.append(item)
			data[item] = 'Completed'
		else:	
			tmp_status = check_status(item)
			if tmp_status == 'MAXIMUM':
				'Reached Max optimization cycle'
				#print(item, '\t\tMaximum opt reached. Restart')
				restart.append(item)
				data[item] = 'Restart'
			elif tmp_status == 'SCF failed to converge':
				'SCF Failed to converge'
				#print(item, '\t\tScf Failed to converge')
				if diverged(calc_dir+'/'+item) == False:
					diverge.append(item)
				else:
					scf_converge.append(item)
				data[item] = 'SCF failed to converge'
			else:
				try:	
					job_ids = job_id_list(calc_dir+'/'+item)
					mutual  = set(running_jobs) & set(job_ids)   
	
					if frozen_atoms() == True:
						'checks atoms are not stuck'
						if len(mutual) != 0:
							'Frozen and Running calculations'
							#print(item, '\t\t', 'Must terminate')
							terminate.append(item)
							for i in mutual:
								terminate_id = terminate_id+' '+str(i)
							data[item] = 'Frozen'
						else:
							'Frozen finished calculations'
							frozen.append(item)	
							#print(item, '\t\t', 'Atoms are frozen. Adjust initial position')
							data[item] = 'Frozen'
					elif len(mutual) != 0:
						'Running calculations (not frozen)'
						#print(item, '\t\t', 'Running')
						Running.append(item)	
						data[item] = 'Running'
					else:	
						#print(item, '\t\tNot sure')
						data[item] = 'Failed'
						not_sure.append(item)
		
				except:
					#print(item, '\t\tNot sure')
					data[item] = 'Failed'
					not_sure.append(item)

if len(restart)>0:
	print('******* Calculations require restart:', len(restart), restart)
if len(frozen)>0:
	print('******* Frozen!:', len(frozen), frozen)
if len(not_started)>0:
	print('******* Not started!:', len(not_started), not_started)
print('******* Scf failed to converge!:',len(scf_converge),  scf_converge)
if len(not_sure)>0:
	print('******* Not sure!:', not_sure)
if len(diverge)>0:
	print('******* Diverge:', len(diverge), diverge)
if len(four_MR)>0:
	print('******* 4MR', len(four_MR), four_MR)

#print('*******\nMust terminate:', len(terminate), terminate)
#print('*******\nRunning:', len(Running), Running)
#print('Completed:', len(completed), completed)

'identify available traj files'
traj = {} #traj[item][thoery type][calc type][status]
for item in data:
	i = item.find('f-')
	if item[i+2:] not in traj:
		traj[item[i+2:]] = {}

'prepare dictionaries to collect data'
traj_copy = deepcopy(traj)
theory_level = ['GGA','hGGA']
calc_type    = ['opt','sp']
results      = ['running','failed','complete']
for t in traj_copy:
	traj[t] = {}
	for th in theory_level:
		traj[t][th] = {}
		for ty in calc_type:
			traj[t][th][ty] = {}
			for r in results:
				traj[t][th][ty][r] = 0
			
'collect completed, failed, running info into the traj dict'
traj_copy = deepcopy(traj)
for t in traj_copy:
	for item in data:
		if item[item.find('f-')+2:] == t:
			'identify calculations of the same reference'

			'level of theory'
			if 'omega' in item:
				theory = 'hGGA'
			elif 'B97' in item:
				theory = 'GGA'

			'calc_type'
			if 'opt' in item:
				calc_type = 'opt'
			elif 'sp' in item:
				calc_type = 'sp'

			'completed/running/failed'
			if data[item] == 'Completed':
				traj[t][theory][calc_type]['complete'] += 1

			elif data[item] == 'Running':
				traj[t][theory][calc_type]['running'] += 1
			else:	
				traj[t][theory][calc_type]['failed']  += 1
	
os.chdir(cwd)
with open("info_status-"+ads+".json", "w") as write_file:
    json.dump(traj, write_file, indent=4)

