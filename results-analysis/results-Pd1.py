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
dir_Pd		    = 'BEA/Pd1'		#name of dir where the calculations are saved
dir_H		    = 'BEA/H'	#name of directory where comensating protons are saved
exchange	    = 'omega'
calc_type	    = 'sp' 		#opt or sp
#traj		    = ['34.traj', '36.traj', '39.traj', '41.traj', '42.traj', '63.traj', '64.traj', '65.traj', '67.traj', '69.traj', '71.traj', '78.traj', '2.traj', '84.traj', '87.traj', '3.traj', '95.traj', '98.traj', '100.traj', '101.traj', '4.traj', '5.traj', '120.traj', '6.traj', '127.traj', '130.traj', '131.traj', '132.traj', '134.traj', '7.traj', '8.traj', '148.traj', '158.traj', '162.traj', '163.traj', '165.traj', '10.traj', '167.traj', '11.traj', '186.traj', '195.traj', '200.traj', '204.traj', '206.traj', '211.traj', '212.traj', '234.traj']	#save the top 5 candidates in those traj files

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

with open(H_data+'/T-atoms.json','r') as read:
	data_T = json.load(read)

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
Al_distance = {}
for ref in references:

	#print(ref)
	#print(len(data_original[ref]['qm_region'])+2,'\n')

	'loop over each reference'
	x_pos, E, first_item, O_d, label = [],[], True, [],[]	#For plotting purposes

	#if os.path.isdir(results_dir+ref) == False:
	#	'create folder if it does not exist'
	#	os.system('mkdir '+results_dir+ref)
	#os.chdir(results_dir+ref)

	for item in references[ref]:
		'each item under reference'
		index = item[0:-5] #remves .traj from the name
		data_output_entry = calc_index(index, data_output, exchange, calc_type) #check corresponding name in data_output
		if data_output_entry != 'none':
			'check calcuation dir is available'
			try:
				if data_output[data_output_entry]['status'] == 'complete':
					'check calc is completed, then copy traj files to new folder'
					#os.system('cp '+calc_dir+'/'+data_output_entry+'/qm-initial.traj '+item[0:-5]+'.traj')
					x_pos.append(int(index))	#x-asis position
					if first_item == True:
						E_ref = data_output[data_output_entry]['energy']
					E.append( (data_output[data_output_entry]['energy']- E_ref)*27.2114 ) #convert from Hartree to e. values of energies for y-axis
					#print(data_output[data_output_entry]['energy'])
					first_item = False
					label.append(index)

					'Al-Al distance'
					atoms = io.read(calc_dir+data_output_entry+'/qm-initial.traj')
					Al_distance[ref], n_Al = Al_Al(atoms)
					'O-O distance'
					#O_O_distance = O_O(atoms, n_Al)
					#O_d.append( round(O_O_distance,3) )
					'Pd-H distance'	
					#Pd_H_distance = Pd_H(atoms, n_Al)
					#Pd_H_d.append( round(Pd_H_distance,3) )
			except:
				pass

	if sorted_plot == True:
		if len(E) >0: #avoid references not calculated yet
			'bar plot (sorted)'
			new_x, new_E, x_pts = sort(x_pos, E)
			#if ref in traj:
			#	'aggreagte resuts to run hGGA calc'
			#	for k in new_x[0:5]:
			#		candidates.append(k)
			print(ref, len(new_x), new_x[0:6])
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
coloring,shade = {},{}
'''Plot rxn energy [bar plot]'''
for entry in minimum:
	if minimum[entry] != '':
		ref_H = data_original[entry]['reference'] #reference for 16 H calculations
		zeolite_H = min_H(ref_H, H_data, calc_type)
		if zeolite_H != 0:
				'zeolite_H == 0 means that some calc are incomplete'
				data_output_entry = calc_index(str(minimum[entry]), data_output, exchange, calc_type) #check corresponding name in data_output
				E_qmmm = data_output[data_output_entry]['energy']
				E_rxn = rxn_energy(E_qmmm, zeolite_H, 1)
				#print(entry, round(E_rxn,2))
			
				#if first_item == True:
				#	E_ref = deepcopy(E_rxn)
				#	first_item = False
				#E.append(E_rxn - E_ref)
				E.append(E_rxn)
				E_label.append(minimum[entry])
				ref_label[minimum[entry]] = entry
				if data_original[entry]['Al'] == 1:
					coloring[entry] = 'y'
				elif data_original[entry]['Al-Al MR'] == [5,5]:
					coloring[entry] = 'r'
					#print(data_original[entry]['Al-Al MR'])	
				elif data_original[entry]['Al-Al MR'] == [6,6]:
					coloring[entry] = 'b'
					#print(data_original[entry]['Al-Al MR'])
				elif data_original[entry]['Al-Al MR'] == [4,4]:
					#print(data_original[entry]['Al-Al MR'])
					coloring[entry] = 'g'
				else:
					coloring[entry] = 'c'
				if data_original[entry]['N'] == 'NNN':
					shade[entry] = '**'
				else:
					shade[entry] = ''
		#else:
		#	#print('{} H calculations are incomplete'.format(entry))
plt.clf()	#clear plot

'Overall Pd Rxn energy plot'
new_x, new_E, x_pts = sort(E_label, E)
#T_atom   = [1119, 1124, 1129, 1134, 1139, 1144, 1149, 1155, 1158]
#for T in T_atom:
for index, item in enumerate(new_E):
		#if ref_label[new_x[index]] in data_T[str(T)]:
		#	plt.bar(x_pts[index], new_E[index],color=coloring[ref_label[new_x[index]]], edgecolor='k', linewidth=4,align='center', alpha=1)
		#else:
		plt.bar(x_pts[index], new_E[index],color=coloring[ref_label[new_x[index]]], hatch=shade[ref_label[new_x[index]]],align='center', alpha=0.9)

		#plt.bar(x_pts[index], new_E[index],color='b', align='center', alpha=1)
		if plt_ref_label == True:
			plt.text(x_pts[index], min(new_E)-0.1, ref_label[new_x[index]], color='k',rotation = 90, fontsize=12)	
			#print(ref_label[new_x[index]])

plt.ylim([-3.1, -2.5])
plt.xticks(x_pts, new_x, rotation = 90)
plt.ylabel('Energy (eV)')
plt.show()
plt.clf()

'Al-Al Distance'
for index, item in enumerate(new_E):
	ref = ref_label[new_x[index]]
	try:
		plt.plot(Al_distance[ref], new_E[index],'ko')
	except:
		pass

plt.xlim([3, 8])
plt.xlabel('Al-Al distance (A)')
plt.ylabel('Pd Rxn Energy (eV)')
plt.show()

