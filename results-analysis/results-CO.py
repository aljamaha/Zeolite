import json, os
import matplotlib.pyplot as plt
from ase import atoms, io
from copy import deepcopy
from functions import *
import numpy as np

'''
Results Analysis
'''

'Inputs'
sorted_plot	    = True		#if True, bar plots of energies is sorted from lowest to highest
plt_ref_label	    = False		#if True, add label of the reference to the overall plot
plotting_overall    = False		#if True, make an overall plot of all results
dir_CO		    = 'Zeolite-CO'	#name of dir where the calculations are saved
dir_H		    = 'CHA-H-aggregate'	#name of directory where comensating protons are saved
CO_gas		    =  -113.320062347	#energy of CO gas
avg		    = True

'Directory names'
data_dir    = '/home/aljama/'+dir_CO+'/data/'			#dir where json data are saved
calc_dir    = '/home/aljama/'+dir_CO+'/calculations/'		#dir where calculations are done
results_dir = '/home/aljama/'+dir_CO+'/results-analysis/' 	#dir where results are to be saved
H_data      = '/home/aljama/'+dir_H+'/data/'			#dir where data for H adsorptions sites are saved

'Load data from json files'
with open(data_dir+"data_output.json", "r") as read_file:
    data_output = json.load(read_file)
with open(data_dir+"data.json", "r") as read_file:
    data_original = json.load(read_file)
with open(H_data+"data_output.json", "r") as read_file:
    data_output_H = json.load(read_file)

def local_ref(ref):
	'accumulate references for H calculations under a traj reference [ref][ref_H]'
	references_H = [] #references for data
	for item in data_original:
		if data_original[item]['reference'] == ref:
			try:
				if data_original[item]['reference_H'] not in references_H:
					references_H.append(data_original[item]['reference_H'])
			except:
				pass
	return references_H

'dictionary where data are saved'
data_H     = {}		#data[ref][ref_H]
data_CO    = {}		#data[ref][ref_H]
references = {} 	#references for data

'Prepare rferences dictionary'
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
'prepare data_CO dictionary'
for ref in references:
	data_CO[ref] = {}
	for ref_H in local_ref(ref):
		data_CO[ref][ref_H] = {}		

'''Accumulate H data'''
for ref in references:
	first_item = True
	data_H[ref] = {}
	for item in references[ref]:
		index = item[0:-5] #remves .traj from the name
		data_output_entry = calc_index(index, data_output_H) #check corresponding name in data_output
		if data_output_entry != 'none':
			if data_output_H[data_output_entry]['status'] == 'complete':	
				#if first_item == True:
				#	E_ref = data_output_H[data_output_entry]['energy']
				#data_H[ref][item] = (data_output_H[data_output_entry]['energy']- E_ref)*27.2114 ) #convert from Hartree
				#first_item = False
				data_H[ref][item] = (data_output_H[data_output_entry]['energy'])*27.2114 

'''Accumulate CO data'''
for ref in references:
	for item in references[ref]:
		if 'adsorbate' in data_original[item]:
			index = item[0:-5] #remves .traj from the name
			data_output_entry = calc_index(index, data_output) #check corresponding name in data_output
			if data_output_entry != 'none':
				if data_output[data_output_entry]['status'] == 'complete':
					ref_in_CO = data_original[item]['reference_H'] #ref_H
					if ref_in_CO not in data_CO[ref]:
						data_CO[ref][ref_in_CO] = {}
					ref_H_E = data_H[ref][ref_in_CO]	
					data_CO[ref][ref_in_CO][item] = (data_output[data_output_entry]['energy'] -  CO_gas - ref_H_E/27.2114) *27.2114

	if avg == True:	
		'save only avg values'
		for item in data_CO[ref]:
			new_e = []
			for e in data_CO[ref][item]:
				new_e.append(data_CO[ref][item][e])
			try:
				avg_e = np.avg(new_e)
				first = True
				for e in data_CO[ref][item]:
					if first == True:
						data_CO[ref][item][e] = avg_e
						first = False
					else:
						del data_CO[ref][item][e]
			except:
				pass
		

	else:
		'save only minimum values'
		for item in data_CO[ref]:
			new_e = []
			for e in data_CO[ref][item]:
				new_e.append(data_CO[ref][item][e])
			try:
				min_e = min(new_e)
				for e in data_CO[ref][item]:
					if data_CO[ref][item][e] != min_e:
						del data_CO[ref][item][e]
			except:
				pass
	
'Plot CO vs. H'
for ref in references:
	x,y = [],[]
	for ref_H in data_H[ref]:
		try:
			H_BE  = data_H[ref][ref_H]
			CO_BE = list(data_CO[ref][ref_H].values())[0]
			label =  list(data_CO[ref][ref_H].keys())[0][0:-5]
			plt.plot(H_BE, CO_BE, 'bo')
			x.append(H_BE), y.append(CO_BE)
			plt.text(H_BE, CO_BE, label)
		except:
			pass
	try:
		print(ref)
		fit = np.polyfit(x,y,1)
		x_fit = np.linspace(min(x),max(x),10)
		y_fit = fit[0]*x_fit + fit[1]
		plt.plot(x_fit,y_fit,'k-',label = 'y = {}x+{}'.format(round(fit[0],2),round(fit[1],2)))
		plt.xlabel('H BE (eV)')
		plt.ylabel('CO BE (eV)')
		plt.legend()
		plt.show()
	except:
		pass
