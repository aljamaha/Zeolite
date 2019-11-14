import os, sys
from matplotlib import pyplot as plt

'''
Objective: Extract energy data from Q-Chem output
Input: python extract-energy.py <output file name>
'''

output_file = sys.argv[1]	#output file name

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

'Plot'
ax = plt.figure(1)
plt.plot(steps, energy,'o-',markersize=10,linewidth=4.0)
plt.xlabel('Step Number', fontsize=18)
plt.ylabel('Energy', fontsize=18)
plt.show()
