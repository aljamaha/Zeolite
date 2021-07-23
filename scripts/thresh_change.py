import os

os.system("ls > tmp")
folders = [line.rstrip('\n') for line in open('tmp')]
cwd = os.getcwd()

for F in folders:
	os.chdir(cwd+'/'+F)
	f = open('opt.in', 'r')
	lines = f.read().split("\n")
	f.close()
	os.system('mv opt.in opt.in-original')
	g = open('opt.in','w')
	for item in lines:
		if 'THRESH' in item:
			g.write('THRESH 13\n')
		else:
			g.write(item+'\n')
	#os.system('cp ~/job.py .')
	#os.system('python job.py opt')
	#os.system('sbatch -J $PWD opt.job')
