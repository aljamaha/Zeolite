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
plotting    	    = False		#if True, plot results for each reference structure
sorted_plot	    = True		#if True, bar plots of energies is sorted from lowest to highest
plt_ref_label	    = False		#if True, add label of the reference to the overall plot
O_n                 = False		#if True, color code plot based on cation-O distance
dir_Pd		    = 'BEA/Pd1/'		#name of dir where the calculations are saved
dir_H		    = 'BEA/H/'	#name of directory where comensating protons are saved
exchange	    = 'omega'
calc_type	    = 'sp' 		#opt or sp

'Directory names'
data_dir    = '/home/aljama/'+dir_Pd+'/data/'			#dir where json data are saved
calc_dir    = '/home/aljama/'+dir_Pd+'/calculations/'		#dir where calculations are done
results_dir = '/home/aljama/'+dir_Pd+'/results-analysis/' 	#dir where results are to be saved
H_data      = '/home/aljama/'+dir_H+'/data/'			#dir where data for H adsorptions sites are saved
candidates  = []						#run hGGA calc on those

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
Pd_H_d, Pd_H_d_all, minimum = {},[],{}	#saves minimum energy for each reference [minimum['3.traj'] = 22.traj]
Al_distance, n_O,n,oxygen_distances= {},{},{},{}
O_Al,O_Si = {},{}
for ref in references:

	'loop over each reference'
	x_pos, E, first_item, O_d, label, candidates = [],[], True, [],[],[]	#For plotting purposes

	for item in references[ref]:
		'each item under reference'
		index = item[0:-5] #remves .traj from the name
		data_output_entry = calc_index(index, data_output, exchange, calc_type) #check corresponding name in data_output
		if data_output_entry != 'none':
			'check calcuation dir is available'
			try:
				if data_output[data_output_entry]['status'] == 'complete':
					'check calc is completed, then copy traj files to new folder'
					x_pos.append(int(index))	#x-asis position
					if first_item == True:
						E_ref = data_output[data_output_entry]['energy']
					E.append( (data_output[data_output_entry]['energy']- E_ref)*27.2114 ) #convert from Hartree to e. values of energies for y-axis
					first_item = False
					label.append(index)

					'Al-Al distance'
					atoms = io.read(calc_dir+data_output_entry+'/qm-initial.traj')
					Al_distance[ref], n_Al = Al_Al(atoms)

					'# of oxygens next to Al'
					if O_n == True:
						O_Al[item], O_Si[item], oxygen_distances[item] = cation_n_O('Pd',atoms,cutoff=2.51)
						#atoms_tmp = io.read(calc_dir+data_output_entry+'/input.xyz')
						#O_Al[item], O_Si[item], oxygen_distances[item] = cation_n_O('Pd',atoms_tmp,cutoff=2.51)

					'O-O distance'
					#O_O_distance = O_O(atoms, n_Al)
					#O_d.append( round(O_O_distance,3) )

					'Pd-H distance'	
					Pd_H_distance = Pd_H(atoms, n_Al)
					Pd_H_d[item] = Pd_H_distance

			except:
				pass

	if sorted_plot == True:
		if len(E) >0: #avoid references not calculated yet
			'bar plot (sorted)'
			new_x, new_E, x_pts = sort(x_pos, E)
			print(ref, len(new_x), new_x[0:13])
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
		zeolite_H = min_H(ref_H, H_data, calc_type)
		if zeolite_H != 0:
			'zeolite_H == 0 means that some calc are incomplete'
			data_output_entry = calc_index(str(minimum[entry]), data_output, exchange, calc_type) #check corresponding name in data_output
			E_qmmm = data_output[data_output_entry]['energy']
			E_rxn = rxn_energy(E_qmmm, zeolite_H, 1)
			
			E.append(E_rxn)
			E_label.append(minimum[entry])
			ref_label[minimum[entry]] = entry
			if data_original[entry]['Al'] == 1:
				coloring[entry] = 'y'
				edge[entry]     = 'y'
			elif data_original[entry]['Al-Al MR'] == [5,5]:
				coloring[entry] = 'r'
				edge[entry]     = 'r'
			elif data_original[entry]['Al-Al MR'] == [6,6]:
				coloring[entry] = 'b'
				edge[entry]     = 'b'
			elif data_original[entry]['Al-Al MR'] == [4,4]:
				coloring[entry] = 'g'
				edge[entry]     = 'g'
			else:
				coloring[entry] = 'c'
				edge[entry]     = 'c'
			if data_original[entry]['N'] == 'NNN':
				shade[entry] = '**'
			else:
				shade[entry] = ''
			#if data_original[entry]['Al MR']['4']>2:
			#	MR4[entry] = True
					
		#else:
		#	#print('{} H calculations are incomplete'.format(entry))
plt.clf()	#clear plot


'Overall Pd Rxn energy plot'
new_x, new_E, x_pts = sort(E_label, E)

k = 0
ref_list = []
for index, item in enumerate(new_E):
	k+=1
	ref = ref_label[new_x[index]]
	ref_list.append(ref)
	name = str(new_x[index])+'.traj'

	'''
		if ref in data_T[str(T)]:
			plt.bar(x_pts[index], new_E[index],color=coloring[ref_label[new_x[index]]], edgecolor='k', linewidth=4,align='center', alpha=1)	
		else:
			plt.bar(x_pts[index], new_E[index],color=coloring[ref_label[new_x[index]]], linewidth=4,align='center', alpha=1)
	'''

	#print(name, O_Al[name], O_Si[name], oxygen_distances[name])

	try:

		if O_n != True:
			#plt.bar(x_pts[index], new_E[index], color=coloring[ref], edgecolor='k', linewidth=4,align='center', alpha=1)
			plt.bar(x_pts[index], new_E[index], color=coloring[ref], hatch=shade[ref], align='center', alpha=0.9)	
		else:
			#if O_Al[name] + O_Si[name] == 4:
			if O_Al[name] == 4:
				plt.bar(x_pts[index], new_E[index], color='g', hatch=shade[ref], align='center', alpha=0.9)
				#plt.bar(x_pts[index], new_E[index], color='k', hatch=shade[ref], edgecolor= edge[ref], align='center', alpha=0.9)
				#print(new_x[index], n[name])
			#elif O_Al[name] + O_Si[name] == 3:
			elif O_Al[name] == 3:
				plt.bar(x_pts[index], new_E[index], color='r', hatch=shade[ref], align='center', alpha=0.9)
				#plt.bar(x_pts[index], new_E[index], color='c', hatch=shade[ref],  edgecolor= edge[ref], align='center', alpha=0.9)
				#print(new_x[index], n[name])
			#elif O_Al[name] + O_Si[name] == 2:
			elif O_Al[name] == 2:
				plt.bar(x_pts[index], new_E[index], color='b', hatch=shade[ref], align='center', alpha=0.9)
			else:
				print(O_Al[name] + O_Si[name])
				plt.bar(x_pts[index], new_E[index], color='y', hatch=shade[ref], align='center', alpha=0.9)
				#plt.bar(x_pts[index], new_E[index], color='y', hatch=shade[ref], edgecolor= edge[ref], align='center', alpha=0.9)

		try:
			plt.text(x_pts[index], np.max(new_E)+0.12, str(O_Si[name]),  color='k',rotation = 90, fontsize=12)
		except:
			pass

		if plt_ref_label == True:
			plt.text(x_pts[index], min(new_E)-0.1, ref,  color='k',rotation = 90, fontsize=12)
			#plt.text(x_pts[index], min(new_E)+0.1, n_O[ref],  color='k',rotation = 90, fontsize=12)
	except:
		print('Failed', name)
		pass
			
plt.ylim([np.min(new_E)-0.1, np.max(new_E)+0.1])
plt.xticks(x_pts, new_x, rotation = 90)
plt.ylabel('Energy (eV)')
plt.show()
plt.clf()

'Pd-H distance plot'
for index, item in enumerate(E_label):
	#print(E_label[index], E[index], index)
	plt.plot(Pd_H_d[str(item)+'.traj'], E[index],'ko')
plt.show()

print('ref list', ref_list)
print('energies', new_E)

'Al-Al Distance'
for index, item in enumerate(new_E):
	ref = ref_label[new_x[index]]

	#try:
	#	if MR4[ref] == True:
	#		plt.plot(Al_distance[ref], new_E[index],'sk',markersize=9)
	#except:
	#	pass

	try:
		#if ref in inaccessible:
		#	plt.plot(Al_distance[ref], new_E[index],'sk',markersize=9)

		if shade[ref] == '**':
			plt.plot(Al_distance[ref], new_E[index],coloring[ref]+'^',markersize=6)
		else:
			plt.plot(Al_distance[ref], new_E[index],coloring[ref]+'o',markersize=6)
	except:
		pass

plt.xlim([3, 8])
plt.xlabel('Al-Al distance (A)')
plt.ylabel('Pd Rxn Energy (eV)')
plt.show()
