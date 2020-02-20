import os
from matplotlib import pyplot as plt

'''
Objective: Extract energy data from Q-Chem output
Input: python extract-energy.py <output file name>
'''

'find output file'
def output():
	'extracts name of the output files'
	out = []
	os.system('ls > tmp')
	files = [line.rstrip('\n') for line in open('tmp')]
	for f in files:
		if '.out' in f:
			out.append(f)
	os.system('rm tmp')
	
	return out

if len(output()) != 1:
	sys.stderr.write("More than one output files in the directory\n")
	quit()
	
output_file = output()[0]

'read output file contect'
f = open(output_file,'r')
lines = f.read().split("\n")
f.close()

step, energy,steps = 0,[],[]

for line in lines:
	'read energy per step'
	if 'Total energy in the final basis set' in line:
		step += 1
		if step == 1:
			E0 = float(line[-14:])
		energy.append(float(line[-14:])- E0)
		steps.append(step)

'write output energy'
f = open('energy.txt','w')
f.write(str(energy[-1]+E0))
f.close()

'Plot'
ax = plt.figure(1)
plt.plot(steps, energy,'o-',markersize=10,linewidth=4.0)
plt.xlabel('Step Number', fontsize=18)
plt.ylabel('Energy', fontsize=18)
plt.show()
