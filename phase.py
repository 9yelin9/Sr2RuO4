# Magnetic Phase Diagram

import os
import re
import time
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

t1 = -4.000

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

ulist = []
nlist = []

inlist = os.listdir(datadir)
inlist.sort()
for innum in range(len(inlist)):

	outlist = os.listdir(datadir+'/'+inlist[innum])
	outlist = [o for o in outlist if o.find('.txt') == -1]

	if not outlist: continue
	else :
		outlist.sort(key=lambda x : x.split('_')[1], reverse=True)
		outfile = inlist[innum]+'/'+outlist[0]

		m  = re.findall('mt[-]?\d+.\d*', outfile)
		if abs(float(m[0].replace('mt', ''))) > 0.01 :
			u = re.findall('U\d+.\d*', outfile)
			n = re.findall('nt\d+.\d*', outfile)

			ulist.append(-float(u[0].replace('U', ''))/t1)
			nlist.append(float(n[0].replace('nt', '')))

### Make diagram
def phase():
	fig = plt.figure()
	plt.plot(nlist, ulist, '.', label='{}'.format(kind.upper()))

	plt.title('Magnetic Phase Diagram ({}/PM)'.format(kind.upper()))
	plt.xlabel('n')
	plt.ylabel('-U/t1')
	plt.yticks(np.arange(0, 5.1, step=1.25))

	plt.legend()
	plt.savefig('/home/9yelin9/sro/diagram/mpd{}_{}.png'.format(kind, runtime))
	plt.show()

### Run
phase()
