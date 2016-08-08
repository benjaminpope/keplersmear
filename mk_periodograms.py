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
	ap.add_argument('ddir', type=str,  help='Input directory name.')
	args = ap.parse_args()
	ddir = args.ddir

	args = ap.parse_args()
	ddir = args.ddir

	stars = {re.search('%s(.+?)_smear_full.csv' % ddir,star).group(1) for star in glob('%s*_smear_full.csv' % ddir)}
	first = clock()

	f = open('do_periodograms.txt','w')

	for star in stars:
		fname = ddir+star
		starttime = clock()
		f.write('python periodograms.py %s' % fname)

	f.close()

	print_time(clock()-first)

	print 'Finished!'