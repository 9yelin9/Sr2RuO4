# Band Structure & Density of States & Fermi Surface Diagram

import os
import re
import time
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

t = time.time()
tm = time.localtime(t)
runtime = '{}{}{}{}'.format(tm.tm_mon, tm.tm_mday, tm.tm_hour, tm.tm_min)

### Options
parser = argparse.ArgumentParser()
parser.add_argument('-k', '--kind', type=str, default='f', choices=['f', 'a'], help='Choose kind of the magnetization\nf : fm, a : afm')
args = parser.parse_args() 

kind = 'fm' if args.kind == 'f' else 'afm'

### Choose file
datadir = '/home/9yelin9/sro/{}/data'.format(kind)

cnt = 0
inlist = os.listdir(datadir)
inlist.sort()
print('\n ♥  INPUT LIST')
for i in range(len(inlist)):
	print('{:2}. {}'.format(i, inlist[i]))
	cnt +=1
if(cnt < 1):
	print('Nothing in here...\n')
	exit()
innum = int(input('\nEnter index of the directory : '))

cnt = 0
outlist = os.listdir(datadir+'/'+inlist[innum])
outlist = [o for o in outlist if o.find('.txt') == -1]
print('\n ♥  OUTPUT LIST')
for i in range(len(outlist)):
	print('{:2}. {}'.format(i, outlist[i]))
	cnt +=1
if(cnt < 1):
	print('Nothing in here...\n')
	exit()
outnum = int(input('\nEnter index of the directory : '))

p  = re.findall('[a-zA-Z]+', inlist[innum])
pv = re.findall('[-]?[0-9]+[.]+[0-9]+', inlist[innum])
pdic = {}
title = ''
for i in range(1, len(pv)):
	pdic[p[i]] = pv[i]
	title += '{}={} '.format(p[i], pv[i])
title += ' {}\n'.format(kind)

p  = re.findall('[a-zA-Z]+', outlist[outnum])
pv = re.findall('[-]?[0-9]+[.]+[0-9]+', outlist[outnum])
for i in range(1, len(pv)):
	pdic[p[i]] = pv[i]
	title += '{}={} '.format(p[i], pv[i])

fname = title
fname = fname.replace(' ', '')
fname = fname.replace('=', '')
fname = fname.replace('\n', '')

### Open file & Make diagram
def band():
	with open(datadir+'/'+inlist[innum]+'/'+outlist[outnum]+'/band.txt', 'r') as data:
		df = pd.read_csv(data, sep='\t')
		df.dropna(axis=1, inplace=True);
		print(df)

	path = [] 
	aup  = []
	adn  = []
	bup  = []
	bdn  = []
	cup  = []
	cdn  = []
			
	for i in range(len(df)):
		path.append(df.values[i][0])
		aup.append(df.values[i][1])
		adn.append(df.values[i][2])
		bup.append(df.values[i][3])
		bdn.append(df.values[i][4])
		cup.append(df.values[i][5])
		cdn.append(df.values[i][6])

	fig = plt.figure(1)
	plt.plot(path, aup, '.', ms=2, color='lightcoral', label=r'$\alpha\uparrow$') 
	plt.plot(path, adn, '.', ms=2, color='tab:red',    label=r'$\alpha\downarrow$')
	plt.plot(path, bup, '.', ms=2, color='lightgreen', label=r'$\beta\uparrow$') 
	plt.plot(path, bdn, '.', ms=2, color='tab:green',  label=r'$\beta\downarrow$') 
	plt.plot(path, cup, '.', ms=2, color='lightblue',  label=r'$\gamma\uparrow$')  
	plt.plot(path, cdn, '.', ms=2, color='tab:blue',   label=r'$\gamma\downarrow$')

	plt.title(title)
	plt.xlabel('Path')
	plt.ylabel('Energy')
	plt.xticks(np.arange(0, len(df)+1, step=float(pdic['k'])), labels=['$\Gamma$', 'M', 'X', '$\Gamma$'])

	plt.legend()
	plt.savefig('/home/9yelin9/sro/diagram/{}band_{}.png'.format(fname, runtime))
	plt.show()

def dos():
	with open(datadir+'/'+inlist[innum]+'/'+outlist[outnum]+'/band.txt', 'r') as data:
		df = pd.read_csv(data, sep='\t')
		df.dropna(axis=1, inplace=True);
		print(df)

	path = [] 
	aup  = []
	adn  = []
	bup  = []
	bdn  = []
	cup  = []
	cdn  = []
			
	for i in range(len(df)):
		path.append(df.values[i][0])
		aup.append(df.values[i][1])
		adn.append(df.values[i][2])
		bup.append(df.values[i][3])
		bdn.append(df.values[i][4])
		cup.append(df.values[i][5])
		cdn.append(df.values[i][6])

	fig = plt.figure(2)
	plt.hist(aup+adn, bins=200, color='tab:red',   histtype='barstacked', density=True, alpha=0.6, label=r'$\alpha$')
	plt.hist(bup+bdn, bins=200, color='tab:green', histtype='barstacked', density=True, alpha=0.6, label=r'$\beta$' )
	plt.hist(cup+cdn, bins=200, color='tab:blue',  histtype='barstacked', density=True, alpha=0.6, label=r'$\gamma$')

	plt.title(title)
	plt.xlabel('Energy')
	plt.ylabel('DOS')

	plt.legend()
	plt.savefig('/home/9yelin9/sro/diagram/{}dos_{}.png'.format(fname, runtime))
	plt.show()

def surface():
	with open(datadir+'/'+inlist[innum]+'/'+outlist[outnum]+'/surface.txt', 'r') as data:
		df = pd.read_csv(data, sep='\t')
		df.dropna(axis=1, inplace=True);
		print(df)

	aupk1 = []
	aupk2 = []
	adnk1 = []
	adnk2 = []

	bupk1 = []
	bupk2 = []
	bdnk1 = []
	bdnk2 = []

	cupk1 = []
	cupk2 = []
	cdnk1 = []
	cdnk2 = []

	for i in range(len(df)):
		if(df.values[i][2] != 0):
			aupk1.append(df.values[i][0])
			aupk2.append(df.values[i][1])
		if(df.values[i][3] != 0):
			adnk1.append(df.values[i][0])
			adnk2.append(df.values[i][1])
		if(df.values[i][4] != 0):
			bupk1.append(df.values[i][0])
			bupk2.append(df.values[i][1])
		if(df.values[i][5] != 0):
			bdnk1.append(df.values[i][0])
			bdnk2.append(df.values[i][1])
		if(df.values[i][6] != 0):
			cupk1.append(df.values[i][0])
			cupk2.append(df.values[i][1])
		if(df.values[i][7] != 0):
			cdnk1.append(df.values[i][0])
			cdnk2.append(df.values[i][1])

	fig = plt.figure()
	plt.plot(aupk1, aupk2, '.', color='lightcoral', label=r'$\alpha\uparrow$')
	plt.plot(adnk1, adnk2, '.', color='tab:red',    label=r'$\alpha\downarrow$')
	plt.plot(bupk1, bupk2, '.', color='lightgreen', label=r'$\beta\uparrow$')
	plt.plot(bdnk1, bdnk2, '.', color='tab:green',  label=r'$\beta\downarrow$')
	plt.plot(cupk1, cupk2, '.', color='lightblue',  label=r'$\gamma\uparrow$')
	plt.plot(cdnk1, cdnk2, '.', color='tab:blue',   label=r'$\gamma\downarrow$')

	plt.title(title)
	plt.xlabel('Path')
	plt.ylabel('Path')
	plt.xticks([-np.pi, 0, np.pi], labels=['$-\pi$', '0', '$\pi$'])
	plt.yticks([-np.pi, 0, np.pi], labels=['$-\pi$', '0', '$\pi$'])

	plt.legend()
	plt.savefig('/home/9yelin9/sro/diagram/{}surface_{}.png'.format(fname, runtime))
	plt.show()

### Run
band()
dos()
surface()
