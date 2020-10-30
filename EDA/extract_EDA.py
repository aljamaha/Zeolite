import os

'Correct double counting in EDA calculations'

QM_name   = 'QM-svp'
QMMM_name = 'eliminate_linking_O'

def Fragment_E(wd):

	output = {}
	f = open(wd+'/opt.out','r')
	lines = f.readlines()
	f.close()
	for index, line in enumerate(lines):
		if 'Fragment Energies' in line:
			n = index
			break

	output['Pd'] = float(lines[n+1][1:-2])
	output['Z']  = float(lines[n+2][1:-2])

	return output

def EDA_results(wd):	

	output = {}
	f = open(wd+'/opt.out','r')
	lines = f.readlines()
	f.close()
	for index, line in enumerate(lines):
		if 'E_frz' in line:
			output['E_frz'] = float(line[line.find('=')+1:-1])
		if 'E_pol' in line:
			output['E_pol'] = float(line[line.find('=')+1:-1])
		if 'E_vct' in line:
			output['E_vct'] = float(line[line.find('=')+1:-1])
		if 'E_int' in line:
			output['E_int'] = float(line[line.find('=')+1:-1])

	return output

cwd = os.getcwd()
QM   = Fragment_E(cwd+'/'+QM_name)
QMMM = Fragment_E(cwd+'/'+QMMM_name)

print('Fragment energies')
print('QM:', QM)
print('QMMMM:', QMMM)

Pd_pol = (QMMM['Pd'] - QM['Pd'])*2625.5 
print('Correction')
print(Pd_pol)

print('EDA - before correction')
EDA_QM   = EDA_results(cwd+'/'+QM_name)
EDA_QMMM = EDA_results(cwd+'/'+QMMM_name)

print('QM: ', EDA_QM,'\n', 'QMMM', EDA_QMMM)
EDA_QMMM['E_frz'] = EDA_QMMM['E_frz'] + Pd_pol
EDA_QMMM['E_int'] = EDA_QMMM['E_int'] + Pd_pol

print('EDA - after correction')
print('QM:', EDA_QM,'\n', 'QMMM', EDA_QMMM)
