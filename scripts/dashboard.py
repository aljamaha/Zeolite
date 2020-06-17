import json, datetime
from matplotlib import pyplot as plt
import numpy as np

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
print('GGA opt completed calculations: {}%'.format(round(total['GGA']['opt']/(GGA_cumulative)*100,1)))
print('GGA sp completed calculations: {}%'.format(round(total['GGA']['sp']/(GGA_cumulative)*100,1)))
print('hGGA opt completed calculations: {}%'.format(round(total['hGGA']['opt']/hGGA_cumulative*100,1)))
print('hGGA sp completed calculations: {}%'.format(round(total['hGGA']['sp']/hGGA_cumulative*100,1)))
print('Total GGA opt so far {} out of {}'.format(total['GGA'], (GGA_cumulative)))

'Day information'
x = datetime.datetime.now()
#summary[x.strftime("%x")] = total	#date of the year
summary[x.strftime("%j")] = total	#day of the year

'dump data'
with open("summary_by_date.json", "w") as write_file:
    json.dump(summary, write_file, indent=4)

d,GGA_opt,GGA_sp,hGGA_opt,hGGA_sp, ref = [],[],[],[],[],[]

'plot data'
ax = plt.figure(1)

for date in summary:
	ref.append(int(date))	#to references all data to first day as 0
	d.append(int(date) - int(ref[0]))
	GGA_opt.append(round(summary[date]['GGA']['opt']/(GGA_cumulative)*100,1))
	GGA_sp.append(summary[date]['GGA']['sp']/(GGA_cumulative)*100)
	hGGA_opt.append(summary[date]['hGGA']['opt']/(hGGA_cumulative)*100)
	hGGA_sp.append(summary[date]['hGGA']['sp']/hGGA_cumulative*100)

print('\n* GGA opt SUMMARY *')
print(d)	#days
print(GGA_opt)	#% of GGA optimized

plt.plot(d, GGA_opt,'bo-',markersize=6,linewidth=3.0,label='GGA opt')
plt.plot([d[0],d[-1]],[(33*32+4*9)/GGA_cumulative*100, (33*32+4*9)/GGA_cumulative*100 ],'k--')
plt.plot(d, GGA_sp,'b^-',markersize=6,linewidth=3.0,label='GGA sp')
plt.plot(d, hGGA_opt,'ro-',markersize=6,linewidth=3.0,label='hGGA opt')
plt.plot(d, hGGA_sp,'r^-',markersize=6,linewidth=3.0,label='hGGA sp')
plt.ylim([0,100])
plt.xlabel('Time (days)', fontsize=16)
plt.ylabel('% of Completed Calculations', fontsize=16)
plt.tick_params(labelsize=13)
plt.legend()
plt.show()
