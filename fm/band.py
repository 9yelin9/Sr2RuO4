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
o1_up = []
o1_dn = []
o2_up = []
o2_dn = []
o3_up = []
o3_dn = []
        
for i in range(len(df)):
	p.append(df.values[i][0])
	o1_up.append(df.values[i][1])
	o1_dn.append(df.values[i][2])
	o2_up.append(df.values[i][3])
	o2_dn.append(df.values[i][4])
	o3_up.append(df.values[i][5])
	o3_dn.append(df.values[i][6])

fig  = plt.figure(figsize=[6, 9])

band = fig.add_subplot(2, 1, 1)
band.plot(p, o1_up, '.', color='lightcoral', label=r'$\alpha\uparrow$')
band.plot(p, o1_dn, '.', color='tab:red',    label=r'$\alpha\downarrow$')
band.plot(p, o2_up, '.', color='lightgreen', label=r'$\beta\uparrow$')
band.plot(p, o2_dn, '.', color='tab:green',  label=r'$\beta\downarrow$')
band.plot(p, o3_up, '.', color='lightblue',  label=r'$\gamma\uparrow$')
band.plot(p, o3_dn, '.', color='tab:blue',   label=r'$\gamma\downarrow$')

band.set_title(title)
band.set_xlabel('Path')
band.set_ylabel('Energy')
band.set_xticks(np.arange(0, len(df)+1, step=float(vv[0])))
band.set_xticklabels(['$\Gamma$', 'M', 'X', '$\Gamma$'])

dos = fig.add_subplot(2, 1, 2)
dos.hist(o1_up, bins=100, histtype='step', color='lightcoral', label=r'$\alpha\uparrow$')  
dos.hist(o1_dn, bins=100, histtype='step', color='tab:red',    label=r'$\alpha\downarrow$')
dos.hist(o2_up, bins=100, histtype='step', color='lightgreen', label=r'$\beta\uparrow$')   
dos.hist(o2_dn, bins=100, histtype='step', color='tab:green',  label=r'$\beta\downarrow$') 
dos.hist(o3_up, bins=100, histtype='step', color='lightblue',  label=r'$\gamma\uparrow$')  
dos.hist(o3_dn, bins=100, histtype='step', color='tab:blue',   label=r'$\gamma\downarrow$')

dos.set_xlabel('Energy')
dos.set_ylabel('Density of states')

plt.legend(loc='lower center', ncol=6, fontsize=9, bbox_to_anchor=(0.5, -0.3))
plt.savefig('/home/9yelin9/mom/fm/diagram/{}'.format(args.name))
plt.show()
