import json, datetime, sys
from matplotlib import pyplot as plt
import numpy as np

try:
	ads = sys.argv[1]
except:
	print('Error! missing adsorbate name')

'load full information on the status of the calculations'
with open('info_status-'+ads+'.json','r') as read:
	traj = json.load(read)

'load summarized data per day'
with open('summary_by_date-'+ads+'.json','r') as read:
	summary = json.load(read)

'total # of calculations needed'
if ads == 'Pd1':
	GGA_cumulative  = 84*32+4*9
	hGGA_cumulative = 84*5+4*9
elif ads == 'H':
	hGGA_cumulative_opt        = 84*16+4*9
	hGGA_cumulative_sp         = 84*5+4*9

'prepare dictionary'
total = {}
total['GGA'], total['hGGA'],total['traj'] = {},{},0
total['GGA']['opt'], total['GGA']['sp'],total['hGGA']['opt'],total['hGGA']['sp'] = 0,0,0,0
completed, failed  = 0,0
hGGA_opt , hGGA_sp, empty = [],[],[]
single_Al = ['1.traj','31.traj','61.traj','93.traj','125.traj','155.traj','185.traj','216.traj','248.traj']	#have a single Al in structure
total['Pd1'], total['H'] = 32,16

'aggreagate data'
for item in traj:
	total['traj'] += 1
	for theory in traj[item]:
		for TYPE in traj[item][theory]:		
			if TYPE == 'opt':
				if traj[item][theory][TYPE]['complete'] > 15:
					completed += 1
			total[theory][TYPE] += traj[item][theory][TYPE]['complete']

'Total failed calculations'
for item in traj:
	for theory in traj[item]:
		for TYPE in traj[item][theory]:
			failed += traj[item][theory][TYPE]['failed']

'print information'
for item in traj:
	if item not in single_Al:
		'single Al are complete'
		if ads == 'Pd1':
			if traj[item]['GGA']['opt']['failed'] + traj[item]['GGA']['opt']['complete'] == total[ads] :
				if traj[item]['GGA']['opt']['complete'] < 16:
					'Too many failed calculations'
					print('WARNING! many failed calc', item, traj[item]['GGA']['opt'])
				if traj[item]['hGGA']['opt']['running'] == 0 and traj[item]['hGGA']['opt']['complete'] == 0 and traj[item]['hGGA']['opt']['failed'] == 0:
					'candidates for hGGA opt'
					hGGA_opt.append(item)
				elif traj[item]['hGGA']['opt']['running'] != 0:
					'hGGA opt is running'
					print('hGGA opt is running', item)
				elif traj[item]['hGGA']['opt']['complete'] > 4:
					'hGGA opt is complete. Evaluate hGGA sp'
					if traj[item]['hGGA']['sp']['running'] == 0 and traj[item]['hGGA']['sp']['complete'] == 0 and traj[item]['hGGA']['sp']['failed'] == 0:
						'Did not start'
						hGGA_sp.append(item)
					elif traj[item]['hGGA']['sp']['running'] != 0:
						'hGGA sp is running'
						empty.append(item)
					elif traj[item]['hGGA']['sp']['complete'] > 4:
						'completed'
						empty.append(item)
					else:
						print('1 or more hGGA sp failed:', item, traj[item]['hGGA']['sp'])
				else:
					print('1 or more hGGA opt failed:', item, traj[item]['hGGA']['opt'])

			else:
				'GGA calc still running ..'
				print('GGA opt still running .. ',item,  traj[item]['GGA']['opt'] )

		elif ads == 'H':		
			if traj[item]['hGGA']['opt']['running'] != 0:
				'hGGA opt is running'
				print('hGGA opt is running', item, traj[item]['hGGA']['opt'])
			elif traj[item]['hGGA']['opt']['failed'] + traj[item]['hGGA']['opt']['complete'] == total[ads] :
				if traj[item]['hGGA']['opt']['complete'] < 12:
					'Too many failed calculations'
					print('WARNING! not all opt completed', item, traj[item]['hGGA']['opt'])
				if traj[item]['hGGA']['sp']['running'] == 0 and traj[item]['hGGA']['sp']['complete'] == 0 and traj[item]['hGGA']['sp']['failed'] == 0:
					'Did not start'
					print('hGGA sp candidate', item, traj[item]['hGGA'])
					hGGA_sp.append(item)
				elif traj[item]['hGGA']['sp']['running'] != 0:
					'hGGA sp is running'
					empty.append(item)
				elif traj[item]['hGGA']['sp']['complete'] > 4:
					'completed'
					empty.append(item)
				else:
					print('1 or more hGGA opt failed:', item, traj[item]['hGGA']['sp'])

			else:
				'hGGA opt still running ..'
				print('hGGA opt still running .. ',item,  traj[item]['hGGA']['opt'] )

if ads == 'Pd1':
	print('\n*hGGA calc:*')
	print('hGGA opt:', hGGA_opt)
	print('hGGA sp:', hGGA_sp)
print('\n*SUMMARY*')
print('Calculations started on {} out of 93 structures'.format(total['traj']))
print('Completed calculations on {} structures'.format(9+completed))
#print('GGA opt completed calculations: {}%'.format(round(total['GGA']['opt']/(GGA_cumulative)*100,1)))
#print('GGA sp completed calculations: {}%'.format(round(total['GGA']['sp']/(GGA_cumulative)*100,1)))
#print('hGGA opt completed calculations: {}%'.format(round(total['hGGA']['opt']/hGGA_cumulative*100,1)))
#print('hGGA sp completed calculations: {}%'.format(round(total['hGGA']['sp']/hGGA_cumulative*100,1)))
#print('Total GGA opt so far {} out of {}'.format(total['GGA'], (GGA_cumulative)))

'Day information'
x = datetime.datetime.now()
#summary[x.strftime("%x")] = total	#date of the year
summary[x.strftime("%j")] = total	#day of the year

'dump data'
with open("summary_by_date-"+ads+".json", "w") as write_file:
    json.dump(summary, write_file, indent=4)

'plot data'
d,GGA_opt,GGA_sp,hGGA_opt,hGGA_sp, ref = [],[],[],[],[],[]

ax = plt.figure(1)

for date in summary:
	ref.append(int(date))	#to references all data to first day as 0
	d.append(int(date) - int(ref[0]))
	if ads == 'Pd1':
		GGA_opt.append(round(summary[date]['GGA']['opt']/(GGA_cumulative)*100,1))
		#GGA_sp.append(summary[date]['GGA']['sp']/(GGA_cumulative)*100)
		hGGA_opt.append(summary[date]['hGGA']['opt']/(hGGA_cumulative)*100)
		hGGA_sp.append(summary[date]['hGGA']['sp']/hGGA_cumulative*100)
	elif ads == 'H':	
		hGGA_opt.append(summary[date]['hGGA']['opt']/(hGGA_cumulative_opt)*100)
		hGGA_sp.append(summary[date]['hGGA']['sp']/hGGA_cumulative_sp*100)

daily_change = []
for index,item in enumerate(GGA_opt):
	if index not in [0,1]:
		daily_change.append(round(GGA_opt[index] - GGA_opt[index-1],1))
try:
	Average = sum(daily_change)/len(daily_change)
	print('Average:      ',round(sum(daily_change)/len(daily_change),2))
except:
	pass

#print('Including failed calc:', round((GGA_opt[-1]*GGA_cumulative/100+failed)/GGA_cumulative*100,1))
if ads == 'Pd1':
	print('% of failed calc:', round(failed/GGA_cumulative*100,1))

ax = plt.figure(1)
GGA_opt = np.array(GGA_opt)
if ads == 'Pd1':
	plt.plot(d, GGA_opt,'bo-',markersize=6,linewidth=3.0,label='GGA opt')
	plt.plot(d[-1], (GGA_opt[-1]*GGA_cumulative/100+failed)/GGA_cumulative*100,'bx',markersize=6,linewidth=3.0)
plt.plot(d, hGGA_opt,'ro-',markersize=6,linewidth=3.0,label='hGGA opt')
plt.plot(d, hGGA_sp,'g^-',markersize=6,linewidth=3.0,label='hGGA sp')
plt.ylim([0,100])
plt.xlabel('Time (days)', fontsize=14)
plt.ylabel('% of Completed Calculations', fontsize=14)
plt.tick_params(labelsize=13)
plt.legend()

ax = plt.figure(2)
try:
	plt.plot([d[2],d[-1]], [Average,Average],'k--')
except:
	pass
plt.bar(d[2:], daily_change,align='center', alpha=0.5)
plt.xlabel('Time (days)', fontsize=16)
plt.ylabel('Daily Change', fontsize=16)
plt.show()
