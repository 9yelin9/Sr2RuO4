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

### Open & Read file
datadir = '/home/9yelin9/mom/fm/data'

flist = os.listdir(datadir)

cnt = 0;
print('\n[ FILE LIST ]')
for i in range(len(flist)):
	if 'surface' in flist[i]:
		print('{:2}. {}'.format(i, flist[i]))
		cnt+=1;

if(cnt < 1):
	print('Nothing in here...\n')
	exit()

fnum = int(input('\nEnter index of the file : '))

with open(datadir+'/{}'.format(flist[fnum]), 'r') as data:
	v  = re.findall('[a-zA-Z]+', flist[fnum])
	vv = re.findall('[-]?[0-9]+.[0-9]+', flist[fnum])

	title = ''
	for i in range(len(vv)):
		title += '{}={} '.format(v[i], vv[i])

	df = pd.read_csv(data, sep='\t')
	df.dropna(axis=1, inplace=True)
	print(df)

### Make Fermi surface diagram
o1up_k1 = []
o1up_k2 = []
o1dn_k1 = []
o1dn_k2 = []

o2up_k1 = []
o2up_k2 = []
o2dn_k1 = []
o2dn_k2 = []

o3up_k1 = []
o3up_k2 = []
o3dn_k1 = []
o3dn_k2 = []

for i in range(len(df)):
	if(df.values[i][2] != 0):
		o1up_k1.append(df.values[i][0])
		o1up_k2.append(df.values[i][1])
	if(df.values[i][3] != 0):
		o1dn_k1.append(df.values[i][0])
		o1dn_k2.append(df.values[i][1])
	if(df.values[i][4] != 0):
		o2up_k1.append(df.values[i][0])
		o2up_k2.append(df.values[i][1])
	if(df.values[i][5] != 0):
		o2dn_k1.append(df.values[i][0])
		o2dn_k2.append(df.values[i][1])
	if(df.values[i][6] != 0):
		o3up_k1.append(df.values[i][0])
		o3up_k2.append(df.values[i][1])
	if(df.values[i][7] != 0):
		o3dn_k1.append(df.values[i][0])
		o3dn_k2.append(df.values[i][1])

fig  = plt.figure()

surface = fig.add_subplot()
surface.plot(o1up_k1, o1up_k2, '.', color = 'lightcoral', label=r'$\alpha\uparrow$')
surface.plot(o1dn_k1, o1dn_k2, '.', color = 'tab:red',    label=r'$\alpha\downarrow$')
surface.plot(o2up_k1, o2up_k2, '.', color = 'lightgreen', label=r'$\beta\uparrow$')
surface.plot(o2dn_k1, o2dn_k2, '.', color = 'tab:green',  label=r'$\beta\downarrow$')
surface.plot(o3up_k1, o3up_k2, '.', color = 'lightblue',  label=r'$\gamma\uparrow$')
surface.plot(o3dn_k1, o3dn_k2, '.', color = 'tab:blue',   label=r'$\gamma\downarrow$')

surface.set_title(title)
surface.set_xlabel('Path')
surface.set_ylabel('Path')
surface.set_xticks([-np.pi, 0, np.pi])
surface.set_yticks([-np.pi, 0, np.pi])
surface.set_xticklabels(['$-\pi$', '0', '$\pi$'])
surface.set_yticklabels(['$-\pi$', '0', '$\pi$'])

plt.legend()
plt.savefig('/home/9yelin9/mom/fm/diagram/{}'.format(args.name))
plt.show()
