import os, sys
from molmod import *

'''
Objective:    convert qchem output to traj files
Input format: python strip-qchem-output.py <qchem-output filename> <xyz initial input filename> <#of qm atoms>
Output:       xyz files under the folder traj
'''

'Modified version from Jeroen Van Der Mynsbrugge original script'

try:
	'input arguments'
	atoms   = sys.argv[3]   #total number of atoms
	with open("data_dir.json", "r") as read_file:
    		data = json.load(read_file)
		atoms = data['total_atoms']
except:
	print('dir_data.json is not found')
	exit()
try:
	'read q-chem output file'
	f = open('opt.out', 'r')
	lines = f.read().split("\n")
	f.close()
except:
	print('opt.out file is not found')
	exit()

if os.path.exists('traj-full atoms') == False:
	os.system('mkdir traj-full-atoms')	#data are saved here

file_number = 0				#for numbering traj files

for index, line in enumerate(lines):
	'identify where qm regions are defined'
	if 'Standard Nuclear Orientation' in line:
		if ' ----------------------------------------------------------------' in lines[index+int(atoms)+3]:
			file_number +=1
			g = open('traj-full-atoms/'+str(file_number)+'.xyz','w')
			g.write(atoms+'\n\n')
			for i in range(index,index+int(atoms)+3):
				if i in [index,1+index,2+index]:
					'avoid empty lines'
					continue
				else:
					g.write(lines[i][11:]+'\n')
			g.close()

os.system('cp traj-full-atoms/'+str(file_number)+'.xyz full-atoms.xyz')
