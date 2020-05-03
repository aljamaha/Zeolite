import os
from ase import io

run  = True
traj = False
cores = '1'
start = 20
end   = 2000
calc_type = 'opt'

cwd = os.getcwd()

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
			pass
		
