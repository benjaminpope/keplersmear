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
mpl.rcParams['savefig.dpi']= 200			 #72 


if __name__ == '__main__':
	ap = ArgumentParser(description='Process smear light curves for periodograms')
	ap.add_argument('star', type=str,  help='Input star name.')
	ap.add_argument('ddir', type=str,  help='Input data directory.')

	args = ap.parse_args()
	star = args.star
	ddir = args.ddir

	first = clock()

	### do high frequencies

	lspmin, lspmax = 1./24.,10

	freqs = np.linspace(1./lspmax, 1./lspmin, 18000)*2.*np.pi

	print 'Doing star %s' % star.replace ("_", " ")
	starttime = clock()

	# read in data
	try:
		print 'Using combined'
		data = Table.read('%s%s_smear_combined.csv' % (ddir,star))
	except:
		data = Table.read('%s%s_smear_full.csv' % (ddir,star))
	good = data['QUARTER']!=0 # quarter 0 is no good
	data = data[good]

	sigs = []
	for q in range(17):
		m = data['QUARTER'] == q

		med, sig = medsig(data['FLUX_CORR_8'][m])

		sigs.append(sig/med)

	sigs = np.array(sigs)
	m2, s2 = medsig(sigs)

	for q in range(17):
		if np.abs(m2-sigs[q]) > (2*s2):
			m = data['QUARTER'] != q
			data = data[m]
			print 'Throwing away quarter', q

	time, flux, corr_flux, filt, quarters = data['BJD'], data['FLUX'],\
	 data['FLUX_CORR_8'], data['FILT'], data['QUARTER']

	thesetimes = time
	thesefluxes = corr_flux

	thesetimes = thesetimes[np.isfinite(thesefluxes)] 
	thesefluxes = thesefluxes[np.isfinite(thesefluxes)]

	# do and plot a lomb-scargle periodogram

	lsp = my_lombscargle(thesetimes,thesefluxes,freqs)

	plt.clf()
	plt.plot(freqs/2./np.pi*11.57,np.sqrt(lsp)*1e6,'k')
	# plt.plot(1./(freqs/2./np.pi),(lsp/lsp.max()),'r')

	fmax = (freqs/2./np.pi*11.57)[np.argmax(lsp)]
	print 'Best period',1./fmax
	print 'Best frequency',fmax
	plt.xlabel(r'Frequency ($\mu$Hz)')
	plt.ylabel('Amplitude (ppm)')
	plt.title('Periodogram of %s Oscillations' % star.replace ("_", " "),y=1.02)

	plt.savefig('%s%s_periodogram_low.png' % (ddir,star))
	print 'Saved periodogram to %s%s_periodogram.png' % (ddir,star)

	# now do low frequency bit

	lspmin, lspmax = 24./24.,100

	freqs = np.linspace(1./lspmax, 1./lspmin, 18000)*2.*np.pi

	print 'Doing star %s low frequencies' % star.replace ("_", " ")
	starttime = clock()

	# read in data
	try:
		print 'Using combined'
		data = Table.read('%s%s_smear_combined.csv' % (ddir,star))
	except:
		data = Table.read('%s%s_smear_full.csv' % (ddir,star))
	good = data['QUARTER']!=0 # quarter 0 is no good
	data = data[good]

	sigs = []
	for q in range(17):
		m = data['QUARTER'] == q

		med, sig = medsig(data['FLUX_CORR_4'][m])

		sigs.append(sig/med)

	sigs = np.array(sigs)
	m2, s2 = medsig(sigs)

	for q in range(17):
		if np.abs(m2-sigs[q]) > (2*s2):
			m = data['QUARTER'] != q
			data = data[m]
			print 'Throwing away quarter', q

	time, flux, corr_flux, filt, quarters = data['BJD'], data['FLUX'],\
	 data['FLUX_CORR_4'], data['FILT'], data['QUARTER']

	thesetimes = time
	thesefluxes = corr_flux

	thesetimes = thesetimes[np.isfinite(thesefluxes)] 
	thesefluxes = thesefluxes[np.isfinite(thesefluxes)]

	# do and plot a lomb-scargle periodogram

	lsp = my_lombscargle(thesetimes,thesefluxes,freqs)

	plt.clf()
	plt.plot(freqs/2./np.pi*11.57,np.sqrt(lsp)*1e6,'k')
	# plt.plot(1./(freqs/2./np.pi),(lsp/lsp.max()),'r')

	fmax = (freqs/2./np.pi*11.57)[np.argmax(lsp)]
	print 'Best period',1./fmax
	print 'Best frequency',fmax
	plt.xlabel(r'Frequency ($\mu$Hz)')
	plt.ylabel('Amplitude (ppm)')
	plt.title('Periodogram of %s Oscillations' % star.replace ("_", " "),y=1.02)

	plt.savefig('%s%s_periodogram_low.png' % (ddir,star))
	print 'Saved periodogram to %s%s_periodogram_low.png' % (ddir,star)

	# ### Do folded light curve
	# print_time(clock()-starttime)

	# print 'Doing planet search'

	# # now do a BLS search

	# lc = Table({'SAP_FLUX':corr_flux,
	# 			'SAP_QUALITY':~np.isfinite(corr_flux),
	# 			'GP_FCOR':corr_flux,
	# 			'BJD':time,
	# 			'GP_TIME':filt})

	# period_range = (0.7,100)

	# bls_results, bls_epoch, bls_dur = bls_fold(20000000,lc, whiten='gp', verbose=True,
	# 	savefigs=False,period_range=period_range,campaign='Nominal', 
	# 		 threshold = 100,figdir ='.')

	# folded = fold(time-bls_epoch,bls_results.bper,shift=0.4)*bls_results.bper
	# plt.clf()
	# plt.plot(folded,lc['GP_FCOR']-lc['GP_TIME'],'.k')
	# plt.axvline(0.4*bls_results.bper)
	# plt.xlabel('Phase')
	# plt.ylabel('Flux')
	# plt.title('%s Folded Light Curve' % (star.replace ("_", " ")),y=1.02)
	# plt.savefig('%s%s_folded.png' % (ddir,star))
	# print 'Saved planet search to %s%s_folded.png' % (ddir,star)

	# plt.xlim(0.3*bls_results.bper,0.5*bls_results.bper)
	# plt.savefig('%s%s_folded_zoom.png' % (ddir,star))
	# print 'Saved planet search to %s%s_folded.png' % (ddir,star)	

	# print_time(clock()-first)

	print 'Finished!'