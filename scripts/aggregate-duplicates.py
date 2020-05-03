import os

'Data'
group = {}
group[1]=['1.traj', '31.traj', '185.traj']
group[2]= ['2.traj', '38.traj', '135.traj', '161.traj']
group[3]= ['3.traj', '5.traj', '37.traj', '39.traj', '66.traj', '97.traj']
group[4]= ['4.traj', '40.traj', '67.traj', '102.traj', '190.traj', '222.traj', '223.traj']
group[5]= ['6.traj', '36.traj', '63.traj', '94.traj', '192.traj', '195.traj', '251.traj', '253.traj', '256.traj', '258.traj']
group[6]= ['7.traj', '35.traj']
group[7]= ['8.traj', '41.traj', '131.traj', '163.traj', '194.traj', '197.traj', '250.traj', '254.traj', '255.traj', '259.traj']
group[8]= ['9.traj', '32.traj', '70.traj', '101.traj', '128.traj', '132.traj', '157.traj', '165.traj']
group[9]= ['10.traj', '34.traj', '71.traj', '100.traj']
group[10]= ['11.traj', '33.traj']
group[11]= ['12.traj', '54.traj', '145.traj', '175.traj']
group[12]= ['13.traj', '55.traj', '87.traj', '115.traj', '200.traj', '201.traj', '205.traj', '210.traj', '231.traj', '233.traj', '238.traj', '240.traj', '245.traj', '247.traj', '261.traj', '277.traj']
group[13]= ['14.traj', '56.traj', '73.traj', '81.traj', '107.traj', '112.traj', '198.traj', '206.traj', '230.traj', '236.traj', '239.traj', '246.traj']
group[14]= ['15.traj', '46.traj', '57.traj', '89.traj', '117.traj', '120.traj', '138.traj', '143.traj', '147.traj', '153.traj', '169.traj', '173.traj', '178.traj', '179.traj']
group[15]= ['16.traj', '21.traj', '52.traj', '59.traj', '75.traj', '78.traj', '82.traj', '84.traj', '105.traj', '110.traj', '121.traj', '122.traj', '140.traj', '141.traj', '149.traj', '150.traj', '151.traj', '170.traj', '172.traj', '181.traj', '183.traj', '184.traj']
group[16]= ['17.traj', '20.traj', '24.traj', '26.traj', '29.traj', '43.traj', '47.traj', '51.traj', '58.traj', '60.traj', '74.traj', '76.traj', '83.traj', '86.traj', '88.traj', '90.traj', '92.traj', '104.traj', '108.traj', '114.traj', '116.traj', '118.traj', '119.traj', '123.traj']
group[16]= ['18.traj', '49.traj']
group[18]= ['19.traj', '53.traj', '72.traj', '106.traj', '209.traj', '212.traj', '263.traj', '269.traj', '271.traj', '276.traj']
group[19]= ['22.traj', '50.traj']
group[20]= ['23.traj', '45.traj']
group[21]= ['25.traj', '42.traj', '142.traj', '148.traj']
group[22]= ['27.traj', '28.traj', '48.traj', '91.traj']
group[23]= ['30.traj', '44.traj', '215.traj', '267.traj', '270.traj']
group[24]= ['61.traj', '93.traj', '216.traj']
group[24]= ['62.traj', '95.traj', '130.traj', '160.traj']
group[26]= ['64.traj', '96.traj', '188.traj', '217.traj', '226.traj']
group[27]= ['65.traj', '103.traj']
group[28]= ['68.traj', '98.traj', '129.traj', '162.traj', '186.traj', '191.traj', '218.traj', '221.traj', '225.traj', '228.traj']
group[29]= ['69.traj', '99.traj', '126.traj', '156.traj']
group[30]= ['77.traj', '109.traj']
group[31]= ['79.traj', '111.traj']
group[32]= ['80.traj', '124.traj']
group[33]= ['85.traj', '113.traj', '139.traj', '168.traj']
group[34]= ['125.traj', '155.traj']
group[35]= ['127.traj', '134.traj', '158.traj', '166.traj']
group[36]= ['133.traj', '164.traj']
group[37]= ['136.traj', '159.traj', '189.traj', '220.traj', '224.traj']
group[38]= ['137.traj', '167.traj']
group[39]= ['144.traj', '176.traj', '203.traj', '235.traj', '243.traj']
group[40]= ['146.traj', '180.traj']
group[41]= ['152.traj', '177.traj']
group[42]= ['154.traj', '171.traj', '174.traj', '182.traj']
group[43]= ['187.traj', '193.traj', '219.traj', '227.traj', '249.traj', '260.traj']
group[44]= ['196.traj', '252.traj', '257.traj']
group[45]= ['199.traj', '202.traj', '229.traj', '232.traj', '242.traj', '262.traj']
group[46]= ['204.traj', '207.traj', '208.traj', '211.traj', '234.traj', '237.traj', '241.traj', '244.traj', '264.traj', '268.traj', '272.traj', '275.traj']
group[47]= ['213.traj', '265.traj', '273.traj']
group[48]= ['214.traj', '266.traj', '274.traj']
group[49]= ['248.traj']

'Inputs'
calc_dir = '/home/aljama/BEA/calculations/'
cwd	 = os.getcwd()

'calculations folders'
os.chdir(calc_dir)
os.system("ls > tmp")
calculations = [line.rstrip('\n') for line in open('tmp')]  #this is the best way


if os.path.exists(cwd+'/group_data') == False:
	os.system('mkdir '+cwd+'/group_data')
for i in group:
	if os.path.exists(cwd+'/group_data/'+str(i)) == False:
		os.system('mkdir '+cwd+'/group_data/'+str(i))
	os.chdir(cwd+'/group_data/'+str(i))
	for traj in group[i]:
		for calc in calculations:
			if traj[0:-5]+'-' == calc[0:len(traj[0:-5])]+'-':
				os.system('cp '+calc_dir+'/'+calc+'/qm-initial.traj '+traj)
