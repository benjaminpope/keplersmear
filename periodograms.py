import numpy as np 
import matplotlib.pyplot as plt 
import re
import glob
from time import time as clock
from k2_epd_george import print_time
from k2_bls_search import *
from keplersmear.smear_tools import *

import matplotlib as mpl 
mpl.style.use('seaborn-colorblind')

if __name__ == '__main__':
	ap = ArgumentParser(description='Process smear light curves for periodograms')
	ap.add_argument('ddir', type=str,  help='Input directory name.')
	args = ap.parse_args()
	ddir = args.ddir

	args = ap.parse_args()
	ddir = args.ddir

	stars = {re.search('%s(.+?)_smear_full.csv' % ddir,star).group(1) for star in glob.glob('%s*_smear_full.csv' % ddir)}
	first = clock()

	lspmin, lspmax = 0.5/24.,100

	freqs = np.linspace(1./lspmax, 1./lspmin, 24000)*2.*np.pi

	for star in stars:
		print 'Doing star %s' % star.replace ("_", " ")
		fname = ddir+star
		starttime = clock()

		# read in data
		data = Table.read('%s%s_smear_full.csv' % (ddir,name))
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

		plt.savefig('%s%s_periodogram.png' % (ddir,star))
		print 'Saved periodogram to %s%s_periodogram.png' % (ddir,star)

		print_time(clock()-starttime)

		print 'Doing planet search'

		## now do a BLS search

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
		plt.savefig('%s%s_folded.png' % (ddir,star))
		print 'Saved planet search to %s%s_folded.png' % (ddir,star)

	print_time(clock()-first)

	print 'Finished!'