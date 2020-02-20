import os, json

'''
Objective: check H and O neighbor list and identify if they have more than needed connecting atoms
'''

try:
	'json data'
	with open("dir_data.json", "r") as read_file:
		data = json.load(read_file)
	qm = data['qm_region']
except:
	print('dir_data.joson does not exist')
	exit()

try:
	'input script'
	f = open('opt.in','r')
	lines = f.readlines() 
	f.close()
except:
	print('opt.in file does not exist')
	exit()

molecule = [] 		#list of lines in moleucles section
stat     = 'fail'

for line in lines:
	if line.strip() == '$molecule':
		#only store info in molecule section
		stat = 'pass'
	if stat == 'pass':
		molecule.append(line.strip().replace('',''))

del molecule[-1], molecule[1], molecule[0]	#del items not containing atom information

for index, item in enumerate(molecule):
	'identify O/H atoms with more than needed connecting atoms'
	entry = item.split('\t')
	if item[0] == 'O':
		if entry[-1] != ' 0' or entry[-2] != ' 0':
			if index in qm:
				print(index, 'in qm region')
			else:
				print(index, 'not in qm region')
	elif item[0] == 'H':	
		if entry[-1] != ' 0' or entry[-2] != ' 0' or entry[-3] != ' 0':	
			if index in qm:
				print(index, 'in qm region')
			else:
				print(index, 'not in qm region')
	elif item[0:2] == 'Si' or item[0:2] == 'Al':
		if entry[-1] == ' 0' or entry[-2] == ' 0' or entry[-3] == ' 0' or entry[-4] == ' 0':
			if index in qm:
				print(index, 'in qm region')	
			else:
				print(index, 'not in qm region')
	else:
		print(item[0:2], entry)

			
		
