import json, os
import matplotlib.pyplot as plt
from ase import atoms, io
from copy import deepcopy
from functions import *

'''
Results Analysis
'''

'Inputs'
plotting    	    = False		#if True, plot results for each reference structure
sorted_plot	    = True		#if True, bar plots of energies is sorted from lowest to highest
plt_ref_label	    = True		#if True, add label of the reference to the overall plot
#plotting_overall    = False		#if True, make an overall plot of all results
dir_Pd		    = 'BEA/Pd1'	#name of dir where the calculations are saved
dir_H		    = 'BEA/H'	#name of directory where comensating protons are saved

'Directory names'
data_dir    = '/home/aljama/'+dir_Pd+'/data/'			#dir where json data are saved
calc_dir    = '/home/aljama/'+dir_Pd+'/calculations/'		#dir where calculations are done
results_dir = '/home/aljama/'+dir_Pd+'/results-analysis/' 	#dir where results are to be saved
H_data      = '/home/aljama/'+dir_H+'/data/'			#dir where data for H adsorptions sites are saved

'Load data from json files'
with open(data_dir+"data_output.json", "r") as read_file:
    data_output = json.load(read_file)

with open(data_dir+"data.json", "r") as read_file:
    data_original = json.load(read_file)

'accumulate reference entries (templates from which calculations were created and run)'
references = {} #references for data

for item in data_original:
	'entries in references dictionary'
	if data_original[item]['reference'] not in references:
		references[data_original[item]['reference']] = []

for ref in references:
	'name of folders that share same reference'
	for item in data_original:
		if data_original[item]['reference'] == ref:
			if ref != item:
				references[ref].append(item)

'accumulate traj files'
Pd_H_d, Pd_H_d_all, minimum = [],[],{}	#saves minimum energy for each reference [minimum['3.traj'] = 22.traj]

for ref in references:

	print(ref)

	'loop over each reference'
	x_pos, E, first_item, O_d, label = [],[], True, [],[]	#For plotting purposes

	#if os.path.isdir(results_dir+ref) == False:
	#	'create folder if it does not exist'
	#	os.system('mkdir '+results_dir+ref)
	#os.chdir(results_dir+ref)

	for item in references[ref]:
		'each item under reference'
		index = item[0:-5] #remves .traj from the name
		data_output_entry = calc_index(index, data_output) #check corresponding name in data_output
		if data_output_entry != 'none':
			'check calcuation dir is available'
			if data_output[data_output_entry]['status'] == 'complete':
				'check calc is completed, then copy traj files to new folder'
				#os.system('cp '+calc_dir+'/'+data_output_entry+'/qm-initial.traj '+item[0:-5]+'.traj')
				x_pos.append(int(index))	#x-asis position
				if first_item == True:
					E_ref = data_output[data_output_entry]['energy']
				E.append( (data_output[data_output_entry]['energy']- E_ref)*27.2114 ) #convert from Hartree to e. values of energies for y-axis
				first_item = False
				label.append(index)

				'Al-Al distance'
				#atoms = io.read(calc_dir+data_output_entry+'/qm-initial.traj')
				#Al_distance, n_Al = Al_Al(atoms)
				'O-O distance'
				#O_O_distance = O_O(atoms, n_Al)
				#O_d.append( round(O_O_distance,3) )
				'Pd-H distance'	
				#Pd_H_distance = Pd_H(atoms, n_Al)
				#Pd_H_d.append( round(Pd_H_distance,3) )


	if sorted_plot == True:
		if len(E) >0: #avoid references not calculated yet
			'bar plot (sorted)'
			new_x, new_E, x_pts = sort(x_pos, E)
			plt.bar(x_pts, new_E, align='center', alpha=1)
			plt.xticks(x_pts, new_x)
			plt.ylabel('Energy (eV)')
			if plotting == True:
				plt.show()

			'save first structure as the minimum energy'
			try:
				minimum[ref] = new_x[0]
			except:
				minimum[ref] = ''
	else:
			'bar plot (not sorted)'
			plt.bar(x_pos, E, align='center', alpha=1)
			plt.xticks(x_pos, x_pos)
			plt.ylabel('Energy (eV)')
			plt.show()

	'''
		plt.plot(O_d, E, 'o', markersize=6)
		plt.xlabel('O-O Distance (A)', fontsize = 10)
		plt.ylabel('Energy (eV)', fontsize = 10)
		for i, lab in enumerate(label):
			plt.text(O_d[i], E[i], lab)
		plt.show()
		plt.plot(Pd_H_d, E, 'o', markersize=6)
		plt.xlabel('Pd-H Distance (A)', fontsize = 10)
		plt.ylabel('Energy (eV)', fontsize = 10)
		for i, lab in enumerate(label):
			plt.text(Pd_H_d[i], E[i], lab)
		plt.show()
	'''

E,E_label,ref_label, first_item = [],[],{}, True

'''Plot rxn energy [bar plot]'''
for entry in minimum:
	if minimum[entry] != '':
		ref_H = data_original[entry]['reference'] #reference for 16 H calculations
		zeolite_H = min_H(ref_H, H_data)
		data_output_entry = calc_index(str(minimum[entry]), data_output) #check corresponding name in data_output
		E_qmmm = data_output[data_output_entry]['energy']
		E_rxn = rxn_energy(E_qmmm, zeolite_H, 1)
		if first_item == True:
			E_ref = deepcopy(E_rxn)
			first_item = False
		E.append(E_rxn - E_ref)
		E_label.append(minimum[entry])
		ref_label[minimum[entry]] = entry
plt.clf()	#clear plot

new_x, new_E, x_pts = sort(E_label, E)
for index, item in enumerate(new_E):
	plt.bar(x_pts[index], new_E[index],color='b', align='center', alpha=1)
	if plt_ref_label == True:
		plt.text(x_pts[index], min(new_E), ref_label[new_x[index]], rotation = 90, fontsize=12)	
plt.xticks(x_pts, new_x, rotation = 90)
plt.ylabel('Energy (eV)')
plt.show()
