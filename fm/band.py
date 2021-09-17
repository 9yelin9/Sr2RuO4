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
p     = [] 
a_up = []
a_dn = []
b_up = []
b_dn = []
c_up = []
c_dn = []
        
for i in range(len(df)):
	p.append(df.values[i][0])
	a_up.append(df.values[i][1])
	a_dn.append(df.values[i][2])
	b_up.append(df.values[i][3])
	b_dn.append(df.values[i][4])
	c_up.append(df.values[i][5])
	c_dn.append(df.values[i][6])

fig  = plt.figure()
band = fig.add_subplot()
band.plot(p, a_up, '.', ms=2, color = 'lightcoral', label=r'$\alpha\uparrow$') 
band.plot(p, a_dn, '.', ms=2, color = 'tab:red',    label=r'$\alpha\downarrow$')
band.plot(p, b_up, '.', ms=2, color = 'lightgreen', label=r'$\beta\uparrow$') 
band.plot(p, b_dn, '.', ms=2, color = 'tab:green',  label=r'$\beta\downarrow$') 
band.plot(p, c_up, '.', ms=2, color = 'lightblue',  label=r'$\gamma\uparrow$')  
band.plot(p, c_dn, '.', ms=2, color = 'tab:blue',   label=r'$\gamma\downarrow$')

band.set_title(title)
band.set_xlabel('Path')
band.set_ylabel('Energy')
band.set_xticks(np.arange(0, len(df)+1, step=float(vv[0])))
band.set_xticklabels(['$\Gamma$', 'M', 'X', '$\Gamma$'])

plt.legend()
plt.savefig('/home/9yelin9/mom/fm/diagram/{}'.format(args.name))
plt.show()
