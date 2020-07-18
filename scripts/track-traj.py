import json, datetime, sys
from matplotlib import pyplot as plt
import numpy as np

data, traj = {},[]
keys = ['H','Pd1']
single_Al = ['1.traj','31.traj','61.traj','93.traj','125.traj','155.traj','185.traj','216.traj','248.traj']	#have a single Al in structure

for key in keys:
	'load full information on the status of the calculations'
	with open('info_status-'+key+'.json','r') as read:
		data[key] = json.load(read)

'Accumulate traj names'
for key in keys:
	for name in data[key]:
		if name not in traj:
			traj.append(name)

'remove single Al'
for item in single_Al:
	traj.remove(item)

for name in traj:
	tmp = {}
	for key in keys:

		'assume it did not start'
		tmp[key] = 'not started'

		if key == 'H':
			try:
				for item in data[key][name]['hGGA']['opt']:
					if data[key][name]['hGGA']['opt'][item] != 0:
						'opt is running ..'
						tmp[key] = 'opt running ..'
			except:
				pass

			if tmp[key] == 'opt running ..':
				if data[key][name]['hGGA']['opt']['complete'] > 12:
					'opt is done'
					tmp[key] = 'opt complete'

			if tmp[key] == 'opt complete':
				for item in data[key][name]['hGGA']['sp']:
					if data[key][name]['hGGA']['sp']['complete'] != 0:
						'sp is running'
						tmp[key] = 'sp is running ..'

			if tmp[key] == 'sp is running':
				if data[key][name]['hGGA']['sp'] > 4:
					'sp is done'
					tmp[key] = 'sp complete'
				

		elif key == 'Pd1':

			try:
				for item in data[key][name]['GGA']['opt']:
					if data[key][name]['GGA']['opt'][item] != 0:
						'GGA opt is running ..'
						tmp[key] = 'GGA opt running ..'
			except:
				pass	

			if tmp[key] == 'GGA opt running ..':
				if data[key][name]['GGA']['opt']['complete'] > 15:
					'GGA opt is done'
					tmp[key] = 'GGA opt complete'

			if tmp[key] == 'GGA opt complete':
				for item in data[key][name]['hGGA']['opt']:
					if data[key][name]['hGGA']['opt']['complete'] != 0:
						'hGGA opt is running'
						tmp[key] = 'hGGA opt running ..'

			if tmp[key] == 'hGGA opt running ..':
				if data[key][name]['hGGA']['opt']['complete'] > 4:
					'hGGA opt is done'
					tmp[key] = 'hGGA opt complete'

			if tmp[key] == 'hGGA opt complete':	
				for item in data[key][name]['hGGA']['sp']:
					if data[key][name]['hGGA']['sp'][item] != 0:
						'hGGA sp running .. '
						tmp[key] = 'hGGA sp running ..'

			if tmp[key] == 'hGGA sp running ..':	
				if data[key][name]['hGGA']['sp']['complete'] > 4:
					'hGGA sp complete'
					tmp[key] = 'hGGA sp done'

			#print(tmp[key], data[key][name]['hGGA']['opt'], name)
			
	print(name, tmp)
