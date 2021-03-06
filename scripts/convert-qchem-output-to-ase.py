import os, sys
from molmod import *

'''
Objective:    convert qchem output to traj files
Input format: python strip-qchem-output.py <qchem-output filename> <xyz initial input filename> <#of qm atoms>
Output:       xyz files under the folder traj
'''

'Modified version from Jeroen Van Der Mynsbrugge original script'

if len(sys.argv) !=2:
	'raise an error when less than two arguments are provided'
	sys.stderr.write("Wrong input format\nInput format: python convert-qchem-output-to-aes.py <#of qm atoms>\n")
	quit()

'input arguments'
qmatoms = sys.argv[1]   #qm atoms (includes also H)
qout    = 'opt.out'

'read q-chem output file'
try:
	f = open(qout, 'r')
	lines = f.read().split("\n")
	f.close()
except:
	print('opt.out file not found')
	exit()

if os.path.exists('traj'):
	d = 'do nothing'
else:
	os.system('mkdir traj')	#data are saved here

file_number = 0		#for numbering traj files

for index, line in enumerate(lines):
	'identify where qm regions are defined'
	if 'Standard Nuclear Orientation' in line:
		if ' ----------------------------------------------------------------' in lines[index+int(qmatoms)+3]:
			file_number +=1
			g = open('traj/'+str(file_number)+'.xyz','w')
			g.write(qmatoms+'\n\n')
			for i in range(index,index+int(qmatoms)+3):
				if i in [index,1+index,2+index]:
					'avoid empty lines'
					continue
				else:
					g.write(lines[i][11:]+'\n')
			g.close()

os.system('cp traj/{}.xyz qm-final.xyz'.format(file_number))
