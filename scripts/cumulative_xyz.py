from ase import io
import os

'accumulate .xyz of multiple qm runs into a single output'

'inputs'
details  = '-opt-omegab97x-d-def2-svp'		#details of the calculations
dir_name = 'output'				#desired name of where xyz files are printed
calc_dir = '/home/aljama/CHA/calculations'	#directory where calculations are done
calc = []					#index of calculations to be written in xyz
for i in range(32,48):
	calc.append(str(i))

'general inputs (dont change)'
number = 0
os.system('mkdir '+dir_name)

'write traj files'
for i in calc:
	number += 1
	os.system('cp '+calc_dir+'/'+i+details+'/qm-final.xyz '+dir_name+'/'+str(number+1)+'.xyz')
