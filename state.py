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
print('\n ♥  INPUT LIST'.format(datadir))
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
print('\n ♥  [{}] OUTPUT LIST'.format(inlist[innum]))
for i in range(len(outlist)):
	print('{:2}. {}'.format(i, outlist[i]))
	cnt +=1
if(cnt < 1):
	print('Nothing in here...\n')
	exit()
outnum = int(input('\nEnter index of the directory : '))

pdic = {}
for i in range(len(inlist[innum].split('_'))):
	pdic[re.sub('[-]?\d+[.]?\d*', '', inlist[innum].split('_')[i])] = re.sub('[a-zA-Z]', '', inlist[innum].split('_')[i])
for j in range(len(outlist[outnum].split('_'))-1):
	pdic[re.sub('[-]?\d+[.]?\d*', '', outlist[outnum].split('_')[j])] = re.sub('[a-zA-Z]', '', outlist[outnum].split('_')[j])

title = '{}={} {}={} {}={} {}\n{}={} {}={} {}={} {}={} {}={}'.format('U', pdic['U'], 'J', pdic['J'], 'k', pdic['k'], kind, 'nt', pdic['nt'], 'mt', pdic['mt'], 'et', pdic['et'], 'mu', pdic['mu'], 'itr', pdic['itr'])

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

	if(kind == 'fm'):
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

		plt.legend(loc='center right', bbox_to_anchor=(1.2, 0.5))
		plt.tight_layout()
		plt.savefig('/home/9yelin9/sro/diagram/{}band_{}.png'.format(fname, runtime))
		plt.show()
	else:
		path = [] 

		a1up  = []
		a1dn  = []
		b1up  = []
		b1dn  = []
		c1up  = []
		c1dn  = []
				
		a2up  = []
		a2dn  = []
		b2up  = []
		b2dn  = []
		c2up  = []
		c2dn  = []

		for i in range(len(df)):
			path.append(df.values[i][0])

			a1up.append(df.values[i][1])
			a1dn.append(df.values[i][2])
			b1up.append(df.values[i][3])
			b1dn.append(df.values[i][4])
			c1up.append(df.values[i][5])
			c1dn.append(df.values[i][6])
              
			a2up.append(df.values[i][7])
			a2dn.append(df.values[i][8])
			b2up.append(df.values[i][9])
			b2dn.append(df.values[i][10])
			c2up.append(df.values[i][11])
			c2dn.append(df.values[i][12])

		fig = plt.figure(1)
		plt.plot(path, a1up, '.', ms=2, color='lightcoral', label=r'$\alpha1\uparrow$') 
		plt.plot(path, a1dn, '.', ms=2, color='tab:red',    label=r'$\alpha1\downarrow$')
		plt.plot(path, b1up, '.', ms=2, color='lightgreen', label=r'$\beta1\uparrow$') 
		plt.plot(path, b1dn, '.', ms=2, color='tab:green',  label=r'$\beta1\downarrow$') 
		plt.plot(path, c1up, '.', ms=2, color='lightblue',  label=r'$\gamma1\uparrow$')  
		plt.plot(path, c1dn, '.', ms=2, color='tab:blue',   label=r'$\gamma1\downarrow$')
                         
		plt.plot(path, a2up, '.', ms=2, color='lightcoral', label=r'$\alpha2\uparrow$') 
		plt.plot(path, a2dn, '.', ms=2, color='darkred',    label=r'$\alpha2\downarrow$')
		plt.plot(path, b2up, '.', ms=2, color='lightgreen', label=r'$\beta2\uparrow$') 
		plt.plot(path, b2dn, '.', ms=2, color='darkgreen',  label=r'$\beta2\downarrow$') 
		plt.plot(path, c2up, '.', ms=2, color='lightblue',  label=r'$\gamma2\uparrow$')  
		plt.plot(path, c2dn, '.', ms=2, color='darkblue',   label=r'$\gamma2\downarrow$')

		plt.title(title)
		plt.xlabel('Path')
		plt.ylabel('Energy')
		plt.xticks(np.arange(0, len(df)+1, step=float(pdic['k'])), labels=['$\Gamma$', 'M', 'X', '$\Gamma$'])

		plt.legend(loc='center right', bbox_to_anchor=(1.25, 0.5))
		plt.tight_layout()
		plt.savefig('/home/9yelin9/sro/diagram/{}band_{}.png'.format(fname, runtime))
		plt.show()

def dos():
	with open(datadir+'/'+inlist[innum]+'/'+outlist[outnum]+'/band.txt', 'r') as data:
		df = pd.read_csv(data, sep='\t')
		df.dropna(axis=1, inplace=True);
		print(df)

	if(kind == 'fm'):
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

		plt.legend(loc='center right', bbox_to_anchor=(1.2, 0.5))
		plt.tight_layout()
		plt.savefig('/home/9yelin9/sro/diagram/{}dos_{}.png'.format(fname, runtime))
		plt.show()
	else:
		path = [] 

		a1up  = []
		a1dn  = []
		b1up  = []
		b1dn  = []
		c1up  = []
		c1dn  = []
		     	
		a2up  = []
		a2dn  = []
		b2up  = []
		b2dn  = []
		c2up  = []
		c2dn  = []
				
		for i in range(len(df)):
			path.append(df.values[i][0])

			a1up.append(df.values[i][1])
			a1dn.append(df.values[i][2])
			b1up.append(df.values[i][3])
			b1dn.append(df.values[i][4])
			c1up.append(df.values[i][5])
			c1dn.append(df.values[i][6])
              
			a2up.append(df.values[i][7])
			a2dn.append(df.values[i][8])
			b2up.append(df.values[i][9])
			b2dn.append(df.values[i][10])
			c2up.append(df.values[i][11])
			c2dn.append(df.values[i][12])

		fig = plt.figure(2)
		plt.hist(a1up+a1dn+a2up+a2dn, bins=200, color='tab:red',    histtype='barstacked', density=True, alpha=0.6, label=r'$\alpha$')
		plt.hist(b1up+b1dn+b2up+b2dn, bins=200, color='tab:green',  histtype='barstacked', density=True, alpha=0.6, label=r'$\beta$' )
		plt.hist(c1up+c1dn+c2up+c2dn, bins=200, color='tab:blue',   histtype='barstacked', density=True, alpha=0.6, label=r'$\gamma$')

		plt.title(title)
		plt.xlabel('Energy')
		plt.ylabel('DOS')

		plt.legend(loc='center right', bbox_to_anchor=(1.25, 0.5))
		plt.tight_layout()
		plt.savefig('/home/9yelin9/sro/diagram/{}dos_{}.png'.format(fname, runtime))
		plt.show()

def surface():
	with open(datadir+'/'+inlist[innum]+'/'+outlist[outnum]+'/surface.txt', 'r') as data:
		df = pd.read_csv(data, sep='\t')
		df.dropna(axis=1, inplace=True);
		print(df)

	if(kind == 'fm'):
		aupp1 = []
		aupp2 = []
		adnp1 = []
		adnp2 = []

		bupp1 = []
		bupp2 = []
		bdnp1 = []
		bdnp2 = []

		cupp1 = []
		cupp2 = []
		cdnp1 = []
		cdnp2 = []

		for i in range(len(df)):
			if(df.values[i][2] != 0):
				aupp1.append(df.values[i][0])
				aupp2.append(df.values[i][1])
			if(df.values[i][3] != 0):
				adnp1.append(df.values[i][0])
				adnp2.append(df.values[i][1])
			if(df.values[i][4] != 0):
				bupp1.append(df.values[i][0])
				bupp2.append(df.values[i][1])
			if(df.values[i][5] != 0):
				bdnp1.append(df.values[i][0])
				bdnp2.append(df.values[i][1])
			if(df.values[i][6] != 0):
				cupp1.append(df.values[i][0])
				cupp2.append(df.values[i][1])
			if(df.values[i][7] != 0):
				cdnp1.append(df.values[i][0])
				cdnp2.append(df.values[i][1])

		fig = plt.figure()
		plt.plot(aupp1, aupp2, '.', color='lightcoral', label=r'$\alpha\uparrow$')
		plt.plot(adnp1, adnp2, '.', color='tab:red',    label=r'$\alpha\downarrow$')
		plt.plot(bupp1, bupp2, '.', color='lightgreen', label=r'$\beta\uparrow$')
		plt.plot(bdnp1, bdnp2, '.', color='tab:green',  label=r'$\beta\downarrow$')
		plt.plot(cupp1, cupp2, '.', color='lightblue',  label=r'$\gamma\uparrow$')
		plt.plot(cdnp1, cdnp2, '.', color='tab:blue',   label=r'$\gamma\downarrow$')

		plt.title(title)
		plt.xlabel('Path')
		plt.ylabel('Path')
		plt.xticks([-np.pi, 0, np.pi], labels=['$-\pi$', '0', '$\pi$'])
		plt.yticks([-np.pi, 0, np.pi], labels=['$-\pi$', '0', '$\pi$'])

		plt.legend(loc='center right', bbox_to_anchor=(1.2, 0.5))
		plt.tight_layout()
		plt.savefig('/home/9yelin9/sro/diagram/{}surface_{}.png'.format(fname, runtime))
		plt.show()
	else:
		a1upp1 = []
		a1upp2 = []
		a1dnp1 = []
		a1dnp2 = []

		b1upp1 = []
		b1upp2 = []
		b1dnp1 = []
		b1dnp2 = []

		c1upp1 = []
		c1upp2 = []
		c1dnp1 = []
		c1dnp2 = []

		a2upp1 = []
		a2upp2 = []
		a2dnp1 = []
		a2dnp2 = []

		b2upp1 = []
		b2upp2 = []
		b2dnp1 = []
		b2dnp2 = []

		c2upp1 = []
		c2upp2 = []
		c2dnp1 = []
		c2dnp2 = []

		for i in range(len(df)):
			if(df.values[i][2] != 0):
				a1upp1.append(df.values[i][0])
				a1upp2.append(df.values[i][1])
			if(df.values[i][3] != 0):
				a1dnp1.append(df.values[i][0])
				a1dnp2.append(df.values[i][1])
			if(df.values[i][4] != 0):
				b1upp1.append(df.values[i][0])
				b1upp2.append(df.values[i][1])
			if(df.values[i][5] != 0):
				b1dnp1.append(df.values[i][0])
				b1dnp2.append(df.values[i][1])
			if(df.values[i][6] != 0):
				c1upp1.append(df.values[i][0])
				c1upp2.append(df.values[i][1])
			if(df.values[i][7] != 0):
				c1dnp1.append(df.values[i][0])
				c1dnp2.append(df.values[i][1])

			if(df.values[i][8] != 0):
				a2upp1.append(df.values[i][0])
				a2upp2.append(df.values[i][1])
			if(df.values[i][9] != 0):
				a2dnp1.append(df.values[i][0])
				a2dnp2.append(df.values[i][1])
			if(df.values[i][10] != 0):
				b2upp1.append(df.values[i][0])
				b2upp2.append(df.values[i][1])
			if(df.values[i][11] != 0):
				b2dnp1.append(df.values[i][0])
				b2dnp2.append(df.values[i][1])
			if(df.values[i][12] != 0):
				c2upp1.append(df.values[i][0])
				c2upp2.append(df.values[i][1])
			if(df.values[i][13] != 0):
				c2dnp1.append(df.values[i][0])
				c2dnp2.append(df.values[i][1])

		fig = plt.figure()
		plt.plot(a1upp1, a1upp2, '.', color='lightcoral', label=r'$\alpha1\uparrow$')
		plt.plot(a1dnp1, a1dnp2, '.', color='tab:red',    label=r'$\alpha1\downarrow$')
		plt.plot(b1upp1, b1upp2, '.', color='lightgreen', label=r'$\beta1\uparrow$')
		plt.plot(b1dnp1, b1dnp2, '.', color='tab:green',  label=r'$\beta1\downarrow$')
		plt.plot(c1upp1, c1upp2, '.', color='lightblue',  label=r'$\gamma1\uparrow$')
		plt.plot(c1dnp1, c1dnp2, '.', color='tab:blue',   label=r'$\gamma1\downarrow$')
		plt.plot(a2upp1, a2upp2, '.', color='lightcoral', label=r'$\alpha2\uparrow$')
		plt.plot(a2dnp1, a2dnp2, '.', color='darkred',    label=r'$\alpha2\downarrow$')
		plt.plot(b2upp1, b2upp2, '.', color='lightgreen', label=r'$\beta2\uparrow$')
		plt.plot(b2dnp1, b2dnp2, '.', color='darkgreen',  label=r'$\beta2\downarrow$')
		plt.plot(c2upp1, c2upp2, '.', color='lightblue',  label=r'$\gamma2\uparrow$')
		plt.plot(c2dnp1, c2dnp2, '.', color='darkblue',   label=r'$\gamma2\downarrow$')

		plt.title(title)
		plt.xlabel('Path')
		plt.ylabel('Path')
		plt.xticks([-np.pi, 0, np.pi], labels=['$-\pi$', '0', '$\pi$'])
		plt.yticks([-np.pi, 0, np.pi], labels=['$-\pi$', '0', '$\pi$'])

		plt.legend(loc='center right', bbox_to_anchor=(1.25, 0.5))
		plt.tight_layout()
		plt.savefig('/home/9yelin9/sro/diagram/{}surface_{}.png'.format(fname, runtime))
		plt.show()

### Run
band()
dos()
surface()
