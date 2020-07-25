import os, json

'checks for every opt folder, there is a sp folder'

'Inputs'
local_dir = 'BEA/Pd1/'
exchange  = 'omega'

'General Inputs'
wd  = '/home/aljama/'
aggregate = []

'calculations available'
os.chdir(wd+local_dir+'/calculations')
os.system("ls > tmp")
folders = [line.rstrip('\n') for line in open('tmp')]  #this is the best way

def comp():
	'check opt.out if calculation has completed'
	tmp = False
	try:
		if os.path.exists('opt.out'):
			os.system('tail opt.out > tmp')
			lines =  [line.rstrip('\n') for line in open('tmp')]
			os.system('rm tmp')
			for line in lines:
				if 'Thank you very much for using Q-Chem' in line:
					tmp = True
					break
	except:
		pass
	return tmp

'check sp calc available for each sp'
for calc in folders:
	if 'opt' in calc and exchange in calc:
		index = calc[0:calc.find('-')]	#calc index
		os.chdir(wd+local_dir+'/calculations/'+calc)
		tmp = False
		if comp() == True:
			for item in folders:
				if item[0:calc.find('-')] == index and 'sp' in item and exchange in item:
					'sp calc is available!'
					tmp = True
					#print(calc, 'Done')
					break
			if tmp == False:
				'no sp calc!'
				print(calc, 'sp not available!')
				aggregate.append(index)

print(aggregate)
				
