import numpy as np
import matplotlib.pyplot as plt

import os
from os.path import exists
from glob import glob
import sys
from time import time as clock

import fitsio
import scipy.io as sio
from argparse import ArgumentParser

from k2_run_specs import get_sc_coords
from SuzPyUtils.norm import *
from SuzPyUtils.multiplot import *
from SuzPyUtils.filter import NIF

import astropy
import astropy.wcs
from astropy.coordinates import SkyCoord, ICRS 
from astropy.wcs.utils import skycoord_to_pixel
import astropy.coordinates as coords
import astropy.constants as const
import astropy.units as units
from astropy.table import Table
from scipy.signal import medfilt
from astropy.time import Time,TimeDelta
import jplephem
import de423

from epic_query import my_MASTRADec

from smear_tools import *
from my_kepffi import MASTRADec, sex2dec
from smearcorrection import *

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

import matplotlib as mpl 
mpl.style.use('seaborn-colorblind')

'''----------------------------------------------------
kepler_smear_lc.py - generate a Kepler smear light curve
----------------------------------------------------'''

if __name__ == '__main__':
	ap = ArgumentParser()
	ap.add_argument('--name', type = str, default = 'test', \
					help = 'star name to save')
	ap.add_argument('--cat-file', type = str, default='kepler_inputs.csv',help = 'list of targets')
	ap.add_argument('--smear-type', type = str, default = None, \
					help = 'which smear to use')
	ap.add_argument('--do-plot', action = 'store_true', default = False, \
					help = 'produce plots')
	ap.add_argument('--do-all', action = 'store_true', default = False, \
					help = 'do all quarters')
	ap.add_argument('--just-stitch',action='store_true', default = False, \
					help='just stitch together quarters without recalculating')

	args = ap.parse_args()
	
	smear_name = lambda s: '' if s is None else str(s)

	name = args.name
	do_plot = args.do_plot 
	smear_type = args.smear_type
	do_all = args.do_all
	cat_file = args.cat_file

	out_dir = 'kepler_smear/'
	cbvdir = '/kepler/kepler/FITS/'
	gap_file = '/users/popeb/K2LC/ben/kepler_bad_cads.csv'

	tab = Table.read(cat_file)

	if do_all:
		names = list(tab['Name'][:])
	else:
		names = [name]

	quarters = range(18)

	starttime = clock()

	for name in names:
		if not args.just_stitch:
			for quarter in quarters:
				try:
					lc = do_target(name,quarter,cat_file=cat_file,out_dir=out_dir,
						cbvdir=cbvdir,smear_type=smear_type,do_plot=do_plot,
						gap_file=gap_file)
					lc.write('%s%s_smear_q%d%s.csv' % (out_dir,name,quarter,smear_name(smear_type)))

					print 'Saved to %s%s_smear_q%d%s.csv' % (out_dir,name,quarter,smear_name(smear_type))

				except:
					print '\nFailed on Quarter %d' % quarter

			print_time(clock()-starttime)

		print 'Now doing the full lightcurve!'

		newlc = stitch_lcs(out_dir,name,smear_name,do_plot=do_plot)

		newlc.write('%s%s_smear_full%s.csv' % (out_dir,name,smear_name(smear_type)))
		print 'Saved lightcurve to %s%s_smear_full%s.csv' % (out_dir,name,smear_name(smear_type))

	
	print_time(clock()-starttime)