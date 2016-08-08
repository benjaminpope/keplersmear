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
	ap.add_argument('star', type=str,  help='Input star name.')

	args = ap.parse_args()
	star = args.star


	first = clock()

	lspmin, lspmax = 1./24.,100

	freqs = np.linspace(1./lspmax, 1./lspmin, 18000)*2.*np.pi

	print 'Doing star %s' % star.replace ("_", " ")
	fname = star
	starttime = clock()

	# read in data
	data = Table.read('%s_smear_full.csv' % (star))
	good = data['QUARTER']!=0 # quarter 0 is no good
	data = data[good]

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

	plt.savefig('%s_periodogram.png' % (star))
	print 'Saved periodogram to %s_periodogram.png' % (star)

	print_time(clock()-starttime)

	print 'Doing planet search'

	# now do a BLS search

	lc = Table({'SAP_FLUX':corr_flux,
				'SAP_QUALITY':~np.isfinite(corr_flux),
				'GP_FCOR':corr_flux,
			   'BJD':time,
			   'GP_TIME':filt})

	period_range = (0.7,500)

	bls_results, bls_epoch, bls_dur = bls_fold(20000000,lc, whiten='gp', verbose=True,
		savefigs=False,period_range=period_range,campaign='Nominal', 
			 threshold = 100,figdir ='.')

	folded = fold(time,bls_results.bper)
	plt.clf()
	plt.plot(folded,lc['GP_FCOR'],'.k')
	plt.xlabel('Phase')
	plt.ylabel('Flux')
	plt.title('%s Folded Light Curve' % star.replace ("_", " "),y=1.02)
	plt.savefig('%s_folded.png' % (star))
	print 'Saved planet search to %s_folded.png' % (star)

	print_time(clock()-first)

	print 'Finished!'