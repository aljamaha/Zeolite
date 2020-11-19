import json, os
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from ase import atoms, io
from copy import deepcopy
from functions import *
import numpy as np
from oxygen_stable import *
import matplotlib.font_manager as font_manager

'''
Results Analysis
'''

'Inputs'
plotting    	    = False		#if True, plot results for each reference structure (individual plots)
sorted_plot	    = True		#if True, bar plots of energies is sorted from lowest to highest
plt_ref_label	    = False		#if True, add label of the reference to the overall plot
O_n                 = False		#if True, color code plot based on cation-O distance
dir_Pd		    = 'BEA/Pd2'		#name of dir where the calculations are saved
dir_H		    = 'BEA/H'		#name of directory where comensating protons are saved
exchange	    = 'omega'
calc_type	    = 'sp' 		#opt or sp
renamed		    = False		#renamed from calc. name to new database name
zeolite		    = 'BEA'

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

'Load new names in database'
with open("/home/aljama/rename_data/"+zeolite+'/'+dir_Pd[-3:]+"_rename.json", "r") as read_file:
    rename = json.load(read_file)

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
O_Al,O_Si, completed_traj = {},{},[]

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
						#O_Al[item], O_Si[item], oxygen_distances[item] = cation_n_O('Pd',atoms,cutoff=2.51)
						O_Al[item], O_Si[item], oxygen_distances[item] = cation_n_O('Pd',atoms,cutoff1=3.0, cutoff2=4.5)
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
			print(ref, len(new_x), new_x[0:5])
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


'complete traj list'
full_list = ['2.traj', '3.traj', '4.traj', '5.traj', '6.traj', '7.traj', '8.traj', '9.traj', '10.traj', '11.traj', '12.traj', '14.traj', '17.traj', '42.traj', '18.traj', '20.traj', '21.traj', '22.traj', '25.traj', '26.traj', '27.traj', '28.traj', '43.traj', '36.traj', '39.traj', '41.traj',  '43.traj', '44.traj', '54.traj', '56.traj', '58.traj', '63.traj', '64.traj', '65.traj', '67.traj', '69.traj', '71.traj', '76.traj', '78.traj', '80.traj', '84.traj', '85.traj', '87.traj', '95.traj', '98.traj', '99.traj', '100.traj', '101.traj', '120.traj', '127.traj', '130.traj', '131.traj', '132.traj', '134.traj', '148.traj', '150.traj', '158.traj', '162.traj', '163.traj', '165.traj', '167.traj', '186.traj', '195.traj', '200.traj', '204.traj', '206.traj', '210.traj', '211.traj', '212.traj', '234.traj', '30.traj', '88.traj', '265.traj']

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
			data_output_entry = calc_index(str(minimum[entry]), data_output, exchange, calc_type) #check corresponding name in data_output
			E_qmmm = data_output[data_output_entry]['energy']
			if 'Pd2' in dir_Pd:
				E_rxn = rxn_energy(E_qmmm, zeolite_H, 2)
			elif 'Pd1' or 'CHA-Pd-agg' in dir_Pd:
				E_rxn = rxn_energy(E_qmmm, zeolite_H, 1)
			E.append(E_rxn)
			E_label.append(minimum[entry])
			ref_label[minimum[entry]] = entry
			if data_original[entry]['Al'] == 1:
				coloring[entry] = 'y'
				edge[entry]     = 'y'
			elif data_original[entry]['Al-Al MR'] == [5,5]:
				coloring[entry] = 'b'
				edge[entry]     = 'b'
			elif data_original[entry]['Al-Al MR'] == [6,6]:
				coloring[entry] = 'r'
				edge[entry]     = 'r'
			elif data_original[entry]['Al-Al MR'] == [4,4]:
				coloring[entry] = 'gray'
				edge[entry]     = 'gray'
			else:
				coloring[entry] = 'g'
				edge[entry]     = 'g'
			if data_original[entry]['N'] == 'NNN':
				shade[entry] = 'w//'
			else:
				shade[entry] = ''
			if entry not in completed_traj:
				completed_traj.append(entry)
					
plt.clf()	#clear plot

'References not yet completed:'
not_yet = []
for item in full_list:
	if item not in completed_traj:
		if item not in not_yet:
			not_yet.append(item)

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

		if O_n != True:
			#plt.bar(x_pts[index], new_E[index], color=coloring[ref], edgecolor='k', linewidth=4,align='center', alpha=1)
			plt.bar(x_pts[index], new_E[index], color=coloring[ref],edgecolor='k', hatch=shade[ref], align='center', alpha=0.9)	
		else:
			#if O_Al[name] + O_Si[name] == 4:
			if int(O_Al[name]) == 4:
				plt.bar(x_pts[index], new_E[index], color='c', hatch=shade[ref], align='center', alpha=0.9)
				#plt.bar(x_pts[index], new_E[index], color='k', hatch=shade[ref], edgecolor= edge[ref], align='center', alpha=0.9)
				#print(new_x[index], n[name])
			#elif O_Al[name] + O_Si[name] == 3:	
			elif O_Al[name] == 3:
				plt.bar(x_pts[index], new_E[index], color='r', hatch=shade[ref], align='center', alpha=0.9)
				#plt.bar(x_pts[index], new_E[index], color='c', hatch=shade[ref],  edgecolor= edge[ref], align='center', alpha=0.9)
				#print(new_x[index], n[name])
			#elif O_Al[name] + O_Si[name] == 2:
			#elif O_Al[name] in [2,1]:
			else:
				plt.bar(x_pts[index], new_E[index], color='b', hatch=shade[ref], align='center', alpha=0.9)
				#else:
				#print(O_Al[name] + O_Si[name])
				#plt.bar(x_pts[index], new_E[index], color='y', hatch=shade[ref], align='center', alpha=0.9)
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

'choice of x-axis label'
renamed_new_x = []
if renamed == True:
	for item in new_x:
		for new_name in rename:
			if str(item) == str(new_name[0:new_name.find('-')]):
				if 'Pd2' in dir_Pd:
					new = rename[new_name].replace('Pd+2','Pd$^{+2}$')
				elif 'Pd1' in dir_Pd:
					new = rename[new_name].replace('Pd+','Pd$^+$')
					new = new.replace('H+','H$^+$')
				#renamed_new_x.append(rename[new_name])
				renamed_new_x.append(new)
	plt.xticks(x_pts, renamed_new_x, rotation = 90,fontsize=12, ha='left')
else:
	plt.xticks(x_pts, new_x, rotation = 90, fontsize=12)


plt.xlim(x_pts[0]-0.5, x_pts[-1]+0.5)
print('Total # of structures', len(x_pts))
plt.subplots_adjust(bottom=0.15)
plt.ylabel('E$\mathrm{_{Pd}}$ (eV)', fontsize=18)
#plt.ylabel(r'E$_{Pd}$ (eV)'+r'\mathrm{_{Pd}}', fontsize=18)
plt.tick_params(axis='y', labelsize=14)
#plt.savefig('Pd1.pdf')
#print('calc names',new_x)
#print('calc energies', new_E)
plt.rcParams["font.family"] = "Times New Roman"

font = font_manager.FontProperties(family='serif', size=16)
if 'Pd2' in dir_Pd:
	custom_lines = [Line2D([0], [0], color='gray', lw=4),
                	Line2D([0], [0], color='b', lw=4),
                	Line2D([0], [0], color='r', lw=4),
                	Line2D([0], [0], color='g', lw=4),]
	plt.legend(custom_lines, ['Al pairs in 4MR', 'Al pairs in 5MR', 'Al pairs in 6MR','Al pairs not in same MR'], loc=7, prop=font) #fontsize=16)
elif 'Pd1' in dir_Pd:	
	custom_lines = [Line2D([0], [0], color='gray', lw=4),
                	Line2D([0], [0], color='b', lw=4),
                	Line2D([0], [0], color='r', lw=4),
                	Line2D([0], [0], color='g', lw=4),
                	Line2D([0], [0], color='y', lw=4)]
	plt.legend(custom_lines, ['Al pairs in 4MR', 'Al pairs in 5MR', 'Al pairs in 6MR','Al pairs not in same MR','Isolated Al'], loc=7, prop=font) #fontsize=16)
plt.show()

plt.clf()

'Pd-H distance plot'
for index, item in enumerate(E_label):
	try:
		ref = data_original[str(item)+'.traj']['reference']
		#print(E_label[index], E[index], index)
		plt.plot(Pd_H_d[str(item)+'.traj'], E[index],coloring[ref]+'o')
	except:
		pass

plt.xlabel('Pd-H distance (A)', fontsize=12)
plt.ylabel('Pd Adsorption Energy (eV)', fontsize=12)
#plt.show()

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
#plt.show()

#print('Top 5 aggregated structures:', top5)
