import json

with open('status_info.json','r') as read:
	traj = json.load(read)

cumulative = 84*32+4*9
total = {}
total['GGA'], total['hGGA'],total['traj'] = 0,0,0
completed  = 0

'print information'
for item in traj:
	total['traj'] += 1
	#print('\n** ', item)
	for theory in traj[item]:
		for TYPE in traj[item][theory]:		
			if TYPE == 'opt':
				if traj[item][theory][TYPE]['complete'] > 20:
					completed += 1
				#print(theory, TYPE, traj[item][theory][TYPE])
				total[theory] += traj[item][theory][TYPE]['complete']
				#print(total)

'print information'
for item in traj:
	print('\n** ', item)
	for theory in traj[item]:
		for TYPE in traj[item][theory]:		
			print(theory, TYPE, traj[item][theory][TYPE])
print('\n')

print('Calculations started on {} out of 93 structures'.format(total['traj']))
print('Completed calculations on {} structures'.format(9+completed))
print('Completed calculations: {}%'.format(round(total['GGA']/cumulative*100,2)))
print('Total GGA opt so far {} out of {}'.format(total['GGA'], cumulative))
