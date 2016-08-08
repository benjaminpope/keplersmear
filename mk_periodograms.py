import numpy as np 
import matplotlib.pyplot as plt 
import re
# import glob
from time import time as clock
from k2_epd_george import print_time
from k2_bls_search import *
from keplersmear.smear_tools import *

import matplotlib as mpl 
mpl.style.use('seaborn-colorblind')
mpl.rcParams['savefig.dpi']= 200             #72 


if __name__ == '__main__':
	ap = ArgumentParser(description='Process smear light curves for periodograms')
	'''-----------------------------------
	Try 
	mpisubcmb "a few hours" 88 /users/popeb/keplersmear/do_periodograms.txt
	-----------------------------------'''
	ap.add_argument('ddir', type=str,  help='Input directory name.')
	args = ap.parse_args()
	ddir = args.ddir

	stars = {re.search('%s(.+?)_smear_full.csv' % ddir,star).group(1) for star in glob('%s*_smear_full.csv' % ddir)}
	first = clock()

	f = open('do_periodograms.txt','w')

	for star in stars:
		starttime = clock()
		f.write('ipython --matplotlib=auto /users/popeb/keplersmear/periodograms.py %s %s > %s%s_output.txt\n' 
			% (star, ddir,ddir,star))
	f.close()

	print_time(clock()-first)

	print 'Finished!'