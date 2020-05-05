#!/home/mgcf/software-ws/anaconda/anaconda3/envs/molmod/bin/python

from ase import io, Atom
import os, pickle, json
from copy import deepcopy
from molmod import *
from functions import *
from qm_region import qm_region
from adsorbate import *
from check_duplicates import *
import time

'''
test check duplicates script (how many duplicates are detected)
'''

dir_name= 'BEA/original'
cwd  	= os.getcwd()
struc_dir = cwd+'/../structures_saved'	#dir to store structures
data_dir        = cwd+'/../data'

with open(data_dir+'/data.json','r') as read:
	data = json.load(read)

'''Detecting duplicate structures'''
data = remove_duplicates(data, struc_dir)


