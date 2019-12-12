import os, sys
from molmod import *

'''
Objective:    convert qchem output to traj files
Input format: python strip-qchem-output.py <qchem-output filename> <xyz initial input filename> <#of qm atoms>
Output:       xyz files under the folder traj
'''

'Modified version from Jeroen Van Der Mynsbrugge original script'

if len(sys.argv) !=4:
	'raise an error when less than two arguments are provided'
	sys.stderr.write("Wrong input format\nInput format: python strip-qchem-output.py <qchem-output file name> <xyz initial input> <total #of atoms>\n")
	quit()

'input arguments'
qout    = sys.argv[1]	#qchem output file name
xyzfile = sys.argv[2]	#xyz coordinates original input file name
atoms   = sys.argv[3]   #total number of atoms

'read q-chem output file'
f = open(qout, 'r')
lines = f.read().split("\n")
f.close()

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
