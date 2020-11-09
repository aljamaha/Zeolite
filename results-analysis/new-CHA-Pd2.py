import json, os
import matplotlib.pyplot as plt
from ase import atoms, io
from copy import deepcopy
from functions import *
import numpy as np
from oxygen_stable import *

'''
Results Analysis
'''

'Inputs'
plotting    	    = False		#if True, plot results for each reference structure (individual plots)
sorted_plot	    = True		#if True, bar plots of energies is sorted from lowest to highest
plt_ref_label	    = False		#if True, add label of the reference to the overall plot
dir_Pd		    = 'CHA-Pd2'		#name of dir where the calculations are saved
dir_H		    = 'CHA-Pd2-H'	#name of directory where comensating protons are saved
exchange	    = 'omega'
calc_type	    = 'sp' 		#opt or sp
renamed		    = False		#renamed from calc. name to new database name

'Directory names'
data_dir    = '/home/aljama/'+dir_Pd+'/data/'			#dir where json data are saved
calc_dir    = '/home/aljama/'+dir_Pd+'/calculations/'		#dir where calculations are done
results_dir = '/home/aljama/'+dir_Pd+'/results-analysis/' 	#dir where results are to be saved
H_data      = '/home/aljama/'+dir_H+'/data/'			#dir where data for H adsorptions sites are saved
candidates,top5  = [],[]						#run hGGA calc on those

'Load data from json files'
with open(data_dir+"data_output.json", "r") as read_file:
    data_output = json.load(read_file)
with open(data_dir+"data.json", "r") as read_file:
    data_original = json.load(read_file)


def sort(x_pos, E):
	'''
	Sort energies of structures from highest to lowest
	Inputs:
		x_pos: label (number) of structures
		E    : energy of the structure
	Outputs:
		new_label: sorted name of each label corresponding to new_E
		new_E    : sorted energies
		x_pts	 : for plotting purposes, from 0 to len(E)
	'''

	E_data, new_label, new_E, x_pts = {},[],[],[]
	E_copy = deepcopy(E)

	for index, item in enumerate(E_copy):
		'define a dict with entries being structure name'
		E_data[x_pos[index]] = E_copy[index]
	E_copy.sort() #sort energies from lowest to highest
	for index, E_item in enumerate(E_copy):
		try:
			x_pts.append(index)			#x-axis points in the plot
			x = list(E_data.values()).index(E_item)
			new_label.append(x_pos[x])		#name of the label of each structure
			new_E.append(E_item)			#sorted energies
		except:
			print('failure in ', index)

	return new_label, new_E, x_pts

def min_H(ref, H_data, calc_type):
	'''
	out of the 16 possible configurations of H sites, identify the lowest energy
	Inputs: ref - name of the original zeolite from which the Pd2+ was created (as well as H2+ by definition)
	Output: minimum energy of the zeolite structure with protons (out of the 16 possibilities)
	'''

	'load data of H adsorbed on zeolite'
	#with open(H_data+'/data.json','r') as read_file:
	#	data_H = json.load(read_file)

	with open(H_data+'/data_output.json','r') as read_file:
		data_H_output = json.load(read_file)

	list_ref = []	#save items that share the same reference

	for item in data_H_output:
		'generate a list of items sharing the same zeolite reference (from which H2+ is created)'
		if '-'+ref in item:
			#for item in data_H:
			#if data_H[item]['reference'] == ref:
			#if item != ref:
			list_ref.append(item)

	min_energy = 0	#define a high energy as starting point

	for item in data_H_output:
		'find energy of each item in list of references'
		#if 'traj' in item:
		#'avoids wrong entries such as incomplete or tmp'
		if item in list_ref:
			if calc_type in item:
				try:
					'only extract calc from sp'
					if data_H_output[item]['energy'] < min_energy:
						a = item
						min_energy =  data_H_output[item]['energy']
				except:
					pass

	return min_energy


'''
'Load new names in database'
with open("/home/aljama/rename_data/"+dir_Pd[-3:]+"_rename.json", "r") as read_file:
    rename = json.load(read_file)
'''

'accumulate reference entries (templates from which calculations were created and run)'
references = {} #references for data
for item in data_original:
	'entries in references dictionary'
	if data_original[item]['reference'] not in references:
		references[data_original[item]['reference']] = []

calc_index = {} #complete name of the folder in the calculation
for ref in references:
	'name of folders that share same reference'
	for item in data_output:
		if '-'+ref in item and 'sp' in item:
			if ref != item:
				if item[0:item.find('-')] not in references[ref]:
					references[ref].append(item[0:item.find('-')])
					calc_index[item[0:item.find('-')]] = item

'accumulate traj files'
minimum = {}	#saves minimum energy for each reference [minimum['3.traj'] = 22.traj]

for ref in references:

	'loop over each reference'
	x_pos, E, first_item, O_d, label, candidates = [],[], True, [],[],[]	#For plotting purposes

	for item in references[ref]:

		data_output_entry = calc_index[item]

		try:
			if data_output[data_output_entry]['status'] == 'complete':
				'check calc is completed, then copy traj files to new folder'
				x_pos.append(item)
				if first_item == True:
					E_ref = data_output[data_output_entry]['energy']
				E.append( (data_output[data_output_entry]['energy']- E_ref)*27.2114 ) 
				first_item = False
				label.append(item)
		except:
			pass

	if sorted_plot == True:
		if len(E) >0: #avoid references not calculated yet
			'bar plot (sorted)'
			new_x, new_E, x_pts = sort(x_pos, E)
			for n in new_x[0:5]:
				top5.append(n)
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
			if plotting  == True:
				plt.show()

'''Plot rxn energy [bar plot]'''
E,E_label,ref_label, first_item = [],[],{}, True
coloring,shade,edge,MR4 = {},{},{},{}

for entry in minimum:
	if minimum[entry] != '':
		ref_H = data_original[entry]['reference'] #reference for 16 H calculations
		try:
			zeolite_H = min_H(ref_H, H_data, calc_type)
		except:
			zeolite_H = 0
		if zeolite_H != 0:
			'zeolite_H == 0 means that some calc are incomplete'
			data_output_entry  = calc_index[minimum[entry]]
			E_qmmm = data_output[data_output_entry]['energy']
			if 'Pd2' in dir_Pd:
				E_rxn = rxn_energy(E_qmmm, zeolite_H, 2)
			elif 'Pd1' or 'CHA-Pd-agg' in dir_Pd:
				E_rxn = rxn_energy(E_qmmm, zeolite_H, 1)
			E.append(E_rxn)
			E_label.append(minimum[entry])
			ref_label[minimum[entry]] = entry
		else:
			print('Failed', entry)
		
plt.clf()		

'Overall Pd Rxn energy plot'
new_x, new_E, x_pts = sort(E_label, E)

k = 0
ref_list = []
for index, item in enumerate(new_E):
	k+=1
	ref = ref_label[new_x[index]]
	ref_list.append(ref)
	name = str(new_x[index])+'.traj'

	try:

		plt.bar(x_pts[index], new_E[index], color='k', align='center', alpha=0.9)	
	except:
		print('Failed', name)

plt.ylim([np.min(new_E)-0.1, np.max(new_E)+0.1])

'''
'choice of x-axis label'
renamed_new_x = []
if renamed == True:
	for item in new_x:
		for new_name in rename:
			if str(item) == str(new_name[0:new_name.find('-')]):
				renamed_new_x.append(rename[new_name])
	plt.xticks(x_pts, renamed_new_x, rotation = 90)
else:
	plt.xticks(x_pts, new_x, rotation = 90)
'''

plt.ylabel('Energy (eV)')
print('calc names',new_x)
print('calc energies', new_E)
plt.show()
plt.clf()

