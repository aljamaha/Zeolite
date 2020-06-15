import json, datetime
from matplotlib import pyplot as plt

'load full information on the status of the calculations'
with open('status_info.json','r') as read:
	traj = json.load(read)

'load summarized data per day'
with open('summary_by_date.json','r') as read:
	summary = json.load(read)

'total # of calculations needed'
GGA_cumulative  = 84*32+4*9
hGGA_cumulative = 84*5+4*9

'prepare dictionary'
total = {}
total['GGA'], total['hGGA'],total['traj'] = {},{},0
total['GGA']['opt'], total['GGA']['sp'],total['hGGA']['opt'],total['hGGA']['sp'] = 0,0,0,0
completed  = 0

'aggreagate data'
for item in traj:
	total['traj'] += 1
	for theory in traj[item]:
		for TYPE in traj[item][theory]:		
			if TYPE == 'opt':
				if traj[item][theory][TYPE]['complete'] > 16:
					completed += 1
			total[theory][TYPE] += traj[item][theory][TYPE]['complete']

'print information'
for item in traj:
	print('\n** ', item)
	for theory in traj[item]:
		for TYPE in traj[item][theory]:		
			print(theory, TYPE, traj[item][theory][TYPE])
print('\n*SUMMARY*\n')
print('Calculations started on {} out of 93 structures'.format(total['traj']))
print('Completed calculations on {} structures'.format(9+completed))
print('GGA opt completed calculations: {}%'.format(round(total['GGA']['opt']/(GGA_cumulative-4*9)*100,1)))
print('GGA sp completed calculations: {}%'.format(round(total['GGA']['sp']/(GGA_cumulative-4*9)*100,1)))
print('hGGA opt completed calculations: {}%'.format(round(total['hGGA']['opt']/hGGA_cumulative*100,1)))
print('hGGA sp completed calculations: {}%'.format(round(total['hGGA']['sp']/hGGA_cumulative*100,1)))
print('Total GGA opt so far {} out of {}'.format(total['GGA'], (GGA_cumulative-4*9)))

'Day information'
x = datetime.datetime.now()
summary[x.strftime("%x")] = total

'dump data'
with open("summary_by_date.json", "w") as write_file:
    json.dump(summary, write_file, indent=4)

'plot data'
ax = plt.figure(1)
for date in summary:
	plt.subplot(2,2,1)
	plt.plot(date, summary[date]['GGA']['opt']/(GGA_cumulative-4*9)*100,'bo',markersize=10,linewidth=4.0)
	plt.ylim([0,100])
	plt.title('GGA opt')
	plt.subplot(2,2,2)
	plt.plot(date, summary[date]['GGA']['sp']/(GGA_cumulative-4*8)*100,'ro-',markersize=10,linewidth=4.0)	
	plt.ylim([0,100])
	plt.title('GGA sp')
	plt.subplot(2,2,3)
	plt.plot(date, summary[date]['hGGA']['opt']/(hGGA_cumulative)*100,'go',markersize=10,linewidth=4.0)
	plt.ylim([0,100])
	plt.title('hGGA opt')
	plt.subplot(2,2,4)
	plt.plot(date, summary[date]['hGGA']['sp']/(hGGA_cumulative)*100,'ko-',markersize=10,linewidth=4.0)
	plt.ylim([0,100])
	plt.title('hGGA sp')

#plt.xlabel('Date', fontsize=18)
#plt.ylabel('% of Completed Calculations', fontsize=18)
plt.show()
