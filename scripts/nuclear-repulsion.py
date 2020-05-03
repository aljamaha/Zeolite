import os
from copy import deepcopy

'''Detects duplicates based on the nuclear repulsion energy'''

'Inputs'
cutoff   = 0.01		#cutoff for comparing NRE 
calc_dir = '/home/aljama/BEA/BEA-qm-repulsion/calculations/'
print_traj = True

'General Inputs'
cwd = os.getcwd()
os.chdir(calc_dir)
os.system("ls > tmp")
folders = [line.rstrip('\n') for line in open('tmp')] 
data = {}

def NRE(pwd):
	'retrieves nuclear repulsion energy under the calculation folder'
	os.chdir(pwd)
	f = open('opt.out','r')
	lines = f.read().split("\n")
	for line in lines:	
		'read nuclear repulsion energies' 
		if 'Nuclear Repulsion Energy ' in line:
			energy = float(line[30:-10])*27.211

	return energy

'Retrieve nuclear repulsion energy for each completed calculation'
for i in folders:
	if 'ref' in i :
		try:
			data[i] = NRE(cwd+'/'+i)
		except:
			print('Calculation did not run for '+i)

'Group duplicates based on nuclear repulsion energy'
duplicate = {}
group = 0
duplicate[group] = []
for item in data:
	for item2 in data:
		if item != item2:
			if abs(data[item] - data[item2]) < cutoff:
				duplicate_copy = deepcopy(duplicate)
				tmp = False
				for j in duplicate_copy:
					if item in duplicate[j] or item2 in duplicate[j]:
						tmp = True
						if item2 not in duplicate[j]:
							duplicate[j].append(item2)
						if item not in duplicate[j]:
							duplicate[j].append(item)
						break
				if tmp == False:
					group+=1
					duplicate[group] = [item,item2]

del duplicate[0]
os.chdir(cwd)

'Print Results'
for item in duplicate:
	print(item, duplicate[item])
	if print_traj == True:
		for index, d in enumerate(duplicate[item]):
			if os.path.exists(str(item)) == False:
				os.system('mkdir '+str(item))
			os.system('cp '+str(duplicate[item][index])+'/qm-initial.traj '+str(item)+'/'+str(duplicate[item][index])+'.traj')
