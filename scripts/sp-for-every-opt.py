import os, json

'checks for every opt folder, there is a sp folder'

'Inputs'
local_dir = 'BEA/H/'
exchange  = 'omega'

'General Inputs'
wd  = '/home/aljama/'
aggregate = []

'calculations available'
os.chdir(wd+local_dir+'/calculations')
os.system("ls > tmp")
folders = [line.rstrip('\n') for line in open('tmp')]  #this is the best way

'load data'
with open(wd+'/'+local_dir+"/data/data.json") as read:
	data = json.load(read)
with open(wd+'/'+local_dir+"/data/data_output.json") as read:
	data_output = json.load(read)

'check sp calc available for each sp'
for calc in folders:
	if 'opt' in calc and exchange in calc:
		index = calc[0:calc.find('-')]	#calc index
		if data_output[calc]['status'] == 'complete':
			tmp = False
			for item in folders:
				if item[0:calc.find('-')] == index and 'sp' in item and exchange in item:
					'sp calc is available!'
					tmp = True
					break
			if tmp == False:
				'no sp calc!'
				if data[index+'.traj']['Al'] == 1:
					print(calc, 'sp not available!')
					aggregate.append(index)

print(aggregate)
				
