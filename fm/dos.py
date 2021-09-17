# Band Structure Diagram

import os
import re
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

### Options
parser = argparse.ArgumentParser()

parser.add_argument('-n', '--name', type=str, default='Figure', help='Name figure')

args = parser.parse_args() 

### Choose file
datadir = '/home/9yelin9/mom/fm/data'

flist = os.listdir(datadir)

cnt = 0
print('\n[ FILE LIST ]')
for i in range(len(flist)):
	if 'band' in flist[i]:
		print('{:2}. {}'.format(i, flist[i]))
		cnt+=1

if(cnt < 1):
	print('Nothing in here...\n')
	exit()

fnum = int(input('\nEnter index of the file : '))

### Open & Read file
with open(datadir+'/{}'.format(flist[fnum]), 'r') as data:
	v  = re.findall('[a-zA-Z]+', flist[fnum])
	vv = re.findall('[-]?[0-9]+.[0-9]+', flist[fnum])

	title = ''
	for i in range(len(vv)):
		title += '{}={} '.format(v[i], vv[i])

	df = pd.read_csv(data, sep='\t')
	df.dropna(axis=1, inplace=True);
	print(df)

### Make band structure diagram
a_up = []
a_dn = []
b_up = []
b_dn = []
c_up = []
c_dn = []
        
for i in range(len(df)):
	a_up.append(df.values[i][1])
	a_dn.append(df.values[i][2])
	b_up.append(df.values[i][3])
	b_dn.append(df.values[i][4])
	c_up.append(df.values[i][5])
	c_dn.append(df.values[i][6])

fig = plt.figure()
pdos = fig.add_subplot()
pdos.hist(a_up+a_dn, bins=200, color='tab:red',   histtype='barstacked', density=True, alpha=0.6, label=r'$\alpha$')
pdos.hist(b_up+b_dn, bins=200, color='tab:green', histtype='barstacked', density=True, alpha=0.6, label=r'$\beta$' )
pdos.hist(c_up+c_dn, bins=200, color='tab:blue',  histtype='barstacked', density=True, alpha=0.6, label=r'$\gamma$')

pdos.set_title(title)
pdos.set_xlabel('Energy')
pdos.set_ylabel('Partial DOS')

plt.savefig('/home/9yelin9/mom/fm/diagram/{}'.format(args.name))
plt.legend()
plt.show()
