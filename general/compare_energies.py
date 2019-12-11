import json
import matplotlib.pyplot as plt
import numpy as np

'''compare the energy of a number of calculations to each other'''

'inputs'
data_dir  = '/home/aljama/CHA/data'	#directory where data are saved
data_name = 'data_output.json'		#name of json file where data are storecalc = []
for item in range(32,49):		#choose here the range of calculations of interest
	calc.append(str(item))

'load data'
with open(data_dir+"/"+data_name, "r") as read_file:
    data = json.load(read_file)

def extract_index(name):
	'returns the index of the calculation name'
	index = name.find('-')

	return name[0:index]

energy = {}	#store energy values and the calc name


ref_check = 0	#uses the first calculation as a reference

for item in data:
	index_name = extract_index(item)
	if index_name in calc:
		'only reference those calculations'
		if ref_check == 0:
			'first calc is our reference'
			ref = data[item]['energy']
			ref_check = 1

		energy[index_name] = (data[item]['energy']-ref)*27.2114 #convert from Hartree to eV

print('Energy: ', energy)

'plot'
x_pos, E = [],[]
for item in energy:
	x_pos.append(int(item))	#x-asis position
	E.append(energy[item])	#values of energies for y-axis

plt.bar(x_pos, E, align='center', alpha=1)
plt.xticks(x_pos, x_pos)
plt.ylabel('Energy (eV)')

plt.show()
