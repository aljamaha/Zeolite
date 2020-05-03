
	'''
	if turn_cutoff == 'on':
		'atoms < cutoff to be included'
		data = O_cutoff(Al_atoms, data, traj, atoms, cutoff)

		'identify Si atoms connected to O atoms in qm region (only ones not accounted for yet)'
		for atom in atoms:
			if atom.symbol == 'O' and atom.index in data[traj]['qm_region']:
				for n_O in N_list[atom.index]:
					if n_O not in data[traj]['qm_region']:
						data[traj]['qm_region'].append(n_O)
	'''

