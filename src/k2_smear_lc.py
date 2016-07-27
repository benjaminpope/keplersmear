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

# from K2fov import fov
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
import wget 

import george
from george import kernels

from smear_tools import *

from my_kepffi import MASTRADec, sex2dec
from epic_query import my_MASTRADec

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

import matplotlib as mpl 
mpl.style.use('seaborn-colorblind')

'''----------------------------------------------------
k2_smear_lc.py - generate a Kepler smear light curve
----------------------------------------------------'''

if __name__ == '__main__':
	ap = ArgumentParser()
	ap.add_argument('epic', type = str, help = 'star name to save')
	ap.add_argument('--campaign', type = str, default='C05',help = 'campaign')
	ap.add_argument('--cat-file', type = str, default='C05_smear.csv',help = 'campaign list')
	ap.add_argument('--smear-type', type = str, default = None, \
					help = 'which smear to use')
	ap.add_argument('--do-plot', action = 'store_true', default = False, \
					help = 'produce plots')
	ap.add_argument('--do-all', action = 'store_true', default = False, \
					help = 'produce plots')


	args = ap.parse_args()

	campaign = args.campaign
	epic = args.epic
	do_plot = args.do_plot 
	smear_type = args.smear_type
	do_all = args.do_all
	cat_file = args.cat_file

	out_dir = '%s_smear/' % campaign

	if do_all:
		tab = Table.read(cat_file)
		epics = list(tab['EPIC'][:])
	else:
		epics = [epic]

	print epics 

	starttime = clock()

	for epic in epics:
		epic = str(epic)
		try:
			print '\nExtracting smear light curve for %s, campaign %s' % (epic,campaign) 

			ra, dec, channel, mod, out, col, row = get_pixel_csv(epic,csv_file=cat_file)
			print 'RA %f, Dec %f' % (ra, dec)

			print 'Found star at mod.out %d.%d, channel %d, row %d, column %d' % \
				(mod,out,channel,row,col)

			col -= 12

			# load the smear data 

			print '\nLoading smear data'
			
			fname = get_smear_file(campaign,mod,out)

			smear = load_smear(fname)
			print 'Smear data loaded'

			## now do MJD - bjd correction
			mjd = np.copy(smear['MJD'][np.isfinite(smear['MJD'])])
			bjd = np.copy(smear['MJD'])
			bjd[np.isfinite(smear['MJD'])] = mjd2bjd(mjd,ra,dec)
			t0 = bjd[np.isfinite(bjd)].min()

			print 'Corrected MJD to BJD'

			if do_plot:
				plt.clf()
				plt.plot(smear['smear_flux'].mean(axis=0),'b')
				plt.plot(smear['vsmear_flux'].mean(axis=0),'g')
				plt.axvline(col,color='r')
				plt.xlabel('pix')
				plt.ylabel('Flux (counts)')
				plt.title('Smear Profile, mod.out %d.%d, %s' % (mod,out,campaign))
				plt.savefig('%sprofile_%d%d_%s.png' % (out_dir,mod,out,campaign))
				print 'Saved smear profile to %sprofile_%d%d_%s.png' % (out_dir,mod,out,campaign)

			# find the star 
			print '\nFitting star centroids'

			starposes = get_centroid_series(smear,col-5,col+5,col=smear_type)

			if do_plot:
				plt.clf()
				plt.plot(bjd-t0,starposes-np.median(starposes),'.k')
				plt.xlabel('BJD - %f' % t0)
				plt.ylabel(r'Delta pix')
				plt.title('%s position, %s' % (epic,campaign))
				plt.savefig('%sstarpos_%s_%s.png' % (out_dir,epic,campaign))
				print 'Saved star positions to %sstarpos_%s_%s.png' % (out_dir,epic,campaign)

			# get background flux 
			print '\nExtracting background flux'

			background = get_background(smear,col=smear_type)

			if do_plot:
				plt.clf()
				plt.plot(bjd-t0,background)
				plt.xlabel('BJD - %f' % t0)
				plt.ylabel('Flux (counts/pix)')
				plt.title('Background flux, mod.out %d.%d, %s' % (mod,out,campaign))
				plt.savefig('%sbackground_%d%d_%s.png' % (out_dir,mod,out,campaign))
				print 'Saved background to %sbackground_%d%d_%s.png' % (out_dir,mod,out,campaign)

			# now extract a lightcurve for several different apertures
			print '\nExtracting aperture photometry'

			widths = [2,3,4,5,6]
			tapers = [2,3,3,3,4]

			flux = np.zeros((5,bjd.shape[0]))
			sigs = np.zeros(5)

			for j in range(5):
				flux[j,:] = extract_lc(smear,starposes,background,widths[j],tapers[j],col=smear_type)
				mm, ss = medsig(flux[j,:])
				sigs[j] = ss/mm

			print 'Best aperture:',np.argmin(sigs)

			if do_plot:
				plt.clf()
				for j in range(5):
					plt.plot(bjd-t0,flux[j,:],'.')
				plt.xlabel('BJD -%f' % t0)
				plt.ylabel('Flux (counts)')
				plt.title('%s Light curve, %s' % (epic,campaign))
				plt.savefig('%slcs_%s_%s.png' % (out_dir,epic,campaign))
				print 'Saved light curves to %slcs_%s_%s.png' % (out_dir,epic,campaign)

				plt.clf()
				plt.plot(bjd-t0,flux[np.argmin(sigs),:]/
					(medsig(flux[np.argmin(sigs),:])[0]),'.k')
				plt.xlabel('BJD -%f' % t0)
				plt.ylabel('Relative Flux')
				plt.title('%s Optimal Light curve, %s' % (epic,campaign))
				plt.savefig('%slc_best_%s_%s.png' % (out_dir,epic,campaign))
				print 'Saved best light curve to %slc_best_%s_%s.png' % (out_dir,epic,campaign)


			# query mast for the nearest object

			try:
				nearby = my_MASTRADec(ra,dec)

				nearby_tab = Table.read(nearby)

				x, y, quality = nearby_tab['POS_CORR1'], nearby_tab['POS_CORR2'],\
					nearby_tab['SAP_QUALITY']

				print '\nSet metadata equal to nearby star'
			except:
				print 'Failed to download a nearby light curve'
				try:
					print 'Expanding search!'
					nearby = my_MASTRADec(ra,dec,darcsec=50000)

					nearby_tab = Table.read(nearby)

					x, y, quality = nearby_tab['POS_CORR1'], nearby_tab['POS_CORR2'],\
						nearby_tab['SAP_QUALITY']
					print '\nSet metadata equal to nearby star'
				except:
					print 'Failed totally to find a nearby light curve'
					quality = ~np.isfinite(flux[0,:])
					x, y = np.zeros_like(bjd), np.zeros_like(bjd)

			smooth_test = NIF(flux[np.argmin(sigs),:],101,13)
			med, sig = medsig(flux[np.argmin(sigs),:]-smooth_test)

			try:
				bads = np.abs(flux[np.argmin(sigs),:]-smooth_test) > (10*sig)
				quality[bads] = 1
			except:
				print 'Failed to filter bad cadences'
				print np.shape(flux)
				print x.shape
				print quality.shape
				print np.size(bads)
				print np.sum(bads),np.size(flux[0,:])
			
			x[np.abs(x)>10.] = np.nan
			y[np.abs(y)>10.] = np.nan
	 
			# now save the lightcurve

			lc = Table({'BJD':bjd,
						'TIME':bjd,
						'MJD':smear['MJD'],
						'CAD': smear['cads'],
						'CADENCENO':smear['cads'],
						'FLUX1': flux[0,:],
						'FLUX2': flux[1,:],
						'FLUX3': flux[2,:],
						'FLUX4': flux[3,:],
						'FLUX5': flux[4,:],
						'SAP_FLUX' : flux[np.argmin(sigs),:],
						'SAP_FLUX_ERR' : flux[np.argmin(sigs),:],
						'SAP_QUALITY':quality,
						'POS_CORR1': x,
						'POS_CORR2': y,
						'BACKGROUND':background,
						'STARPOS':starposes
						})

			smear_name = lambda s: '' if s is None else ('_'+str(s))
			lc.write('%s%s%s_smear%s.csv' % (out_dir,campaign,epic,smear_name(smear_type)))
			lc.write('%s%s%s_smear%s.fits' % (out_dir,campaign,epic,smear_name(smear_type)),
				overwrite=True)

			print '\nSaved to %s%s%s_smear%s.csv' % (out_dir,campaign,epic,smear_name(smear_type))
			print '\nSaved to %s%s%s_smear%s.fits too!' % (out_dir,campaign,epic,smear_name(smear_type))
		except:
			print '\nFailed on EPIC %s' % epic

	print_time(clock()-starttime)