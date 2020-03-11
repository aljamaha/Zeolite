from ase import io, Atoms

def qm_structure(data, struc_dir):
	'''
	prints structure of only the qm region
	inputs : data dictionary
	outputs: qm structures
	'''
	for item in data:
		new = Atoms('Ni', [(0, 0, 0)],cell=[1, 1, 1]) #dummy structure
		atoms = io.read(item)
		for j in data[item]['qm_region']:
			new.append(atoms[j])
		del new[0]
		new.write(item[:-5]+'-qm.traj')
