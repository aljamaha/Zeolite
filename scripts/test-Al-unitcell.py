import os
from ase import io

start = 1
end   = 278
cwd = '/home/aljama/BEA/original/calculations'

for i in range(start,end):
	i = str(i)
	try:
			os.chdir(cwd+'/'+i+'-opt-omegab97x-d-def2-svp-ref-'+i+'.traj')	
			atoms = io.read('input.xyz')
			for atom in atoms:
				if atom.symbol == 'Al':
					if atom.index > 192:
						print('Warning!', i)
	except:
			print('Failed!', i)
		
