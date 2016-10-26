import numpy as np
import matplotlib.pyplot as plt

import os
from os.path import exists
from glob import glob
import sys
from time import time as clock
import json
import urllib
import fitsio
import scipy.io as sio
from argparse import ArgumentParser

from k2_run_specs import get_sc_coords
from SuzPyUtils.norm import *
from SuzPyUtils.multiplot import *
from SuzPyUtils.filter import NIF
from scipy.signal import savgol_filter

import astropy
import astropy.wcs
from k2_epd_george import print_time
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
from scipy.signal import lombscargle
from epic_query import my_MASTRADec

import re
import george
from george import kernels

from my_kepffi import sex2dec
# from smearcorrection import *
from kepsys import cbv
from k2sc.cdpp import cdpp


def centroid(thisflux,thispos):
    return np.sum(thisflux*thispos)/np.sum(thisflux)

###----------------------------------------------
###----------------------------------------------


def medsig(array):
    '''
    Return median and outlier-robust estimate of standard deviation
    (1.48 x median of absolute deviations).
    '''
    l = np.isfinite(array)
    if sum(l) == 0:
        return np.nan, np.nan
    if sum(l) == 1:
        return array[l], np.nan
    med = np.median(array[l])
    sig = 1.48 * np.median(abs(array[l] - med))
    return med, sig


###----------------------------------------------
###----------------------------------------------

def sex2dec(ra,dec):
    '''Copied from Tom Barclay's PyKE'''

    ra = re.sub('\s+','|',ra.strip())
    ra = re.sub(':','|',ra.strip())
    ra = re.sub(';','|',ra.strip())
    ra = re.sub(',','|',ra.strip())
    ra = re.sub('-','|',ra.strip())
    ra = ra.split('|')
    outra = (float(ra[0]) + float(ra[1]) / 60 + float(ra[2]) / 3600) * 15.0

    dec = re.sub('\s+','|',dec.strip())
    dec = re.sub(':','|',dec.strip())
    dec = re.sub(';','|',dec.strip())
    dec = re.sub(',','|',dec.strip())
    # dec = re.sub('-.','|',dec.strip())
    dec = dec.split('|')
    if float(dec[0]) > 0.0:
        outdec = float(dec[0]) + float(dec[1]) / 60 + float(dec[2]) / 3600
    else:
        outdec = float(dec[0]) - float(dec[1]) / 60 - float(dec[2]) / 3600

    return outra, outdec


###----------------------------------------------
###----------------------------------------------

def cosbell_filter(width,taper):
    flat = np.ones(width+taper)
    base = np.hanning(taper)
    flat[0:taper/2] = base[0:taper/2]
    flat[-taper/2:] = base[-taper/2:]
    return flat

###----------------------------------------------
###----------------------------------------------

def cosbell(x,pos,width,taper,subsample=4):

    # subsample
    subsampled = np.arange(x.size*subsample)

    # make template

    flat_top = np.zeros_like(subsampled).astype('float64')

    flat_top[(pos*subsample-(subsample/2)*width-(subsample/2)*taper):(pos*subsample+(subsample/2)*width+(subsample/2)*taper)] =\
     cosbell_filter(width*subsample,taper*subsample)

    # rebin
    final = np.zeros_like(x).astype('float64')

    for j in range(subsample):
        final += flat_top[j::subsample]

    return final

###----------------------------------------------
###----------------------------------------------

def supergaussian(x,pos,width):
    return np.exp(-((x-pos)/width)**4.)

###----------------------------------------------
###----------------------------------------------

def supergaussian_bin(x,pos,width,subsample=8):
    subsampled = np.arange(x.size*subsample)
    near = np.abs(subsampled/float(subsample)-pos) < 4*width
    dummy = np.zeros_like(subsampled).astype('float64')

    dummy[near] = np.exp(-((subsampled[near]/float(subsample)-pos)/width)**4.)
    final = np.zeros_like(x).astype('float64')

    for j in range(subsample):
        final += dummy[j::subsample]

    return final

###----------------------------------------------
###----------------------------------------------

def get_num(quarter):
    '''Get the file number for each quarter'''
    if quarter == 0:
        return 2009131105131
    elif quarter == 1:
        return 2009166043257
    elif quarter == 2:
        return 2009259160929
    elif quarter == 3:
        return 2009350155506
    elif quarter == 4: # two sets of collateral - this is released 2014-12-17 16:15:42
        return 2010078095331
    # elif quarter == 4: # two sets of collateral - this is released 2012-04-25 13:41:56
    #   return 2010078170814
    elif quarter == 5:
        return 2010174085026
    elif quarter == 6:
        return 2010265121752
    elif quarter == 7:
        return 2010355172524
    elif quarter == 8:
        return 2011073133259
    elif quarter == 9:
        return 2011177032512
    elif quarter == 10:
        return 2011271113734
    elif quarter == 11:
        return 2012004120508
    elif quarter == 12:
        return 2012088054726
    elif quarter == 13:
        return 2012179063303
    elif quarter == 14:
        return 2012277125453
    elif quarter == 15:
        return 2013011073258
    elif quarter == 16:
        return 2013098041711
    elif quarter == 17:
        return 2013131215648

###----------------------------------------------
###----------------------------------------------

def censor_bad_cads(lc,quarter,gap_file='kepler_bad_cads.csv'):
    bad_cads = Table.read(gap_file)
    try: # not all have bad cadences
        index = np.where(bad_cads['Quarter']==quarter)
        time = lc['BJD']
        start, stop = float(bad_cads['Start'][index]),float(bad_cads['Stop'][index])
        bad = np.where((time>=start) & (time<=stop))
        print '\nCensoring data between',start, stop

        if quarter == 11: # there are two jumps!
            print 'Censoring additional cadences'
            also_bad = np.where((time >= 55926) & (time <=55926.5))
            other_bad = np.where((time >= 55929.3) & (time <=55929.6)) 
            bad = np.concatenate((bad[0],also_bad[0],other_bad[0]))
            # bad = np.concatenate((bad[0],other_bad[0]))

        if quarter == 12: # there are two jumps!
            print 'Censoring additional cadences'
            also_bad = np.where((time >= 55948) & (time <=55952))
            other_bad = np.where((time >= 55953) & (time <=55956)) 
            bad = np.concatenate((bad[0],also_bad[0],other_bad[0]))
            # bad = np.concatenate((bad[0],other_bad[0]))

        for key in lc.keys():
            if 'FLUX' in key:
                lc[key][bad] = np.nan
        print 'Censored %s bad cadences' % np.size(bad)

    except:
        print 'No bad cadences in list'
    return lc

###----------------------------------------------
###----------------------------------------------

def load_kepler(fname,gap_file='kepler_bad_cads.csv'):
    lc = Table.read(fname)
    quarter = int(re.search('q\d{1,2}',fname).group(0)[1:])

    try:
        newlc = censor_bad_cads(lc,quarter,gap_file=gap_file)
        lc = newlc.copy()
    except:
        pass

    return lc

###----------------------------------------------
###----------------------------------------------

def movingaverage(interval, window_size):
    window= np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'same')

###----------------------------------------------
###----------------------------------------------

def fold(time, period, origo=0.0, shift=0.0, normalize=True,  clip_range=None):
    """Folds the given data over a given period.

    from Hannu

    Parameters
    ----------
    
      time  
      period      
      origo 
      shift 
      normalize   
      clip_range  
    Returns
    -------
      phase 
    """
    
    tf = ((time - origo)/period + shift) % 1.

    if not normalize:
        tf *= period
        
    if clip_range is not None:
        mask = np.logical_and(clip_range[0]<tf, tf<clip_range[1])
        tf = tf[mask], mask
    return tf

###----------------------------------------------
###----------------------------------------------

def my_lombscargle(times,fluxes,freqs,norm=True):

    mask = np.isfinite(fluxes)*np.isfinite(times)

    ndata = np.sum(mask)

    thesetimes, thesefluxes = np.copy(times[mask]).byteswap().newbyteorder().astype('float64'),\
     np.copy(fluxes[mask]).byteswap().newbyteorder().astype('float64')

    medflux = np.median(thesefluxes)
    if medflux == 0:
        medflux = 1

    if norm == True:
        thesefluxes -= medflux

    lsp = lombscargle(thesetimes,thesefluxes,freqs)*4./ndata/medflux

    # if norm == True:
    #    lsp /= medflux

    return lsp

###----------------------------------------------
###----------------------------------------------

def my_stft(times,fluxes,freqs,nbins=64):

    ndata = np.shape(fluxes)[0]
    nfreqs = np.shape(freqs)[0]
    nstep = np.ceil(ndata/np.float(nbins))

    stft_amp = np.zeros((nfreqs,nbins))

    for j in range(nbins):
        good = np.isfinite(times[(j*nstep):((j+1)*nstep)])*np.isfinite(fluxes[(j*nstep):((j+1)*nstep)])
        if np.sum(good) < 0.1*nstep:
            continue
        stft_amp[:,j] = my_lombscargle(times[(j*nstep):((j+1)*nstep)],
            fluxes[(j*nstep):((j+1)*nstep)],freqs)

    return stft_amp

###----------------------------------------------
###----------------------------------------------

def k2data_smear(lc,flux_type='SAP_FLUX'):
    time = np.copy(lc['BJD'][:])
    try:
        cadence = np.copy(lc['CAD'][:])
    except:
        cadence = np.arange(np.size(time))
    quality = np.copy(lc['SAP_QUALITY'][:])
    fluxes = np.copy(lc[flux_type][:])
    try:
        errors = np.copy(lc['SAP_FLUX_ERR'][:])
    except:
        errors = 0.001*fluxes
    x, y = np.copy(lc['POS_CORR1'][:]), np.copy(lc['POS_CORR2'][:])
    return time, cadence, quality, fluxes, errors, x, y

###----------------------------------------------
###----------------------------------------------

def detrend_smear(lc,epic=None,flux_type='SAP_FLUX',flare_sigma=5,flare_erosion=5,
    ls_min_period=0.2,ls_max_period=20.,ls_min_power=100.,splits=[57170.],
    de_npop=100,de_niter=150,de_max_time=120.,outlier_sigma=5,outlier_mwidth=5,
    tr_nrandom=300,tr_nblocks=6,tr_bspan=50):
    if epic==None:
        epic = '20000000'

    dataset = K2Data(epic,*k2data_smear(lc,flux_type=flux_type))

    # mask flares 
    # dataset.mask_flares(flare_sigma, flare_erosion)

    # do a Lomb-Scargle

    # dataset.search_for_periodicity(ls_min_period,ls_max_period,ls_min_power)

    if dataset.is_periodic:
        kernel = QuasiPeriodicKernel(period=dataset.ls_period)
        print 'Found periodicity p = {:7.2f} (power {:7.2f} > {:7.2f})'\
            .format(dataset.ls_period, dataset.ls_power, ls_min_power)
    else:
        kernel = BasicKernel()
        print 'No strong periodicity found: using an aperiodic kernel'

    ## initialize a detrender

    tstart = clock()
    inputs = np.transpose([dataset.time,dataset.x,dataset.y])
    detrender = Detrender(dataset.fluxes.ravel(), inputs, mask=np.isfinite(dataset.fluxes.ravel()),
                          splits=splits, kernel=kernel, tr_nrandom=tr_nrandom,
                          tr_nblocks=tr_nblocks, tr_bspan=tr_bspan)
    de = DiffEvol(detrender.neglnposterior, kernel.bounds, de_npop)

    if isinstance(kernel, QuasiPeriodicKernel):
        de._population[:,2] = np.clip(normal(dataset.ls_period, 0.1*dataset.ls_period, size=de.n_pop),
                                      ls_min_period, ls_max_period)

    print 'Starting global hyperparameter optimisation using DE'

    tstart_de = clock()

    for i,r in enumerate(de(de_niter)):
        print 'DE iteration %3i -ln(L) %4.1f' % (i, de.minimum_value)
        tcur_de = clock()
        if ((de._fitness.ptp() < 3) or (tcur_de - tstart_de > de_max_time)) and (i>2):
            break
    print 'DE finished in %i seconds'%  (tcur_de-tstart_de)
    print 'DE minimum found at: %s' %  np.array_str(de.minimum_location, precision=3, max_line_width=250)
    print 'DE -ln(L) %4.1f' % (de.minimum_value)

    print 'Starting local hyperparameter optimisation'

    try:
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=RuntimeWarning, append=True)
            pv, warn = detrender.train(de.minimum_location)
    except ValueError as e:
        print 'Local optimiser failed, %s' % e
        return

    print 'Local minimum found at: %s' % np.array_str(pv, precision=3)

    print 'Masking outliers'
    detrender.mask_outliers(pv=pv)

    print 'Computing time and position trends'
    tr_time, tr_position = detrender.predict(pv, components=True)

    cdpp_r = cdpp(detrender.data.masked_time, detrender.data.masked_flux)
    cdpp_c = cdpp(detrender.data.unmasked_time, detrender.data.unmasked_flux-tr_time-tr_position, exclude=~detrender.data.mask)

    print 'Raw CDPP %6.3f' % cdpp_r
    print 'Detrended CDPP %6.3f' % cdpp_c

    print_time(clock()-tstart)

    lc2 = Table({'SAP_FLUX':detrender.data.unmasked_flux,
                 'SAP_QUALITY':np.ones_like(detrender.data.unmasked_flux),
                 'BJD':detrender.data.unmasked_time,
                 'TIME':detrender.data.unmasked_time,
                 'GP_TIME':tr_time,
                 'GP_POS':tr_position,
                 'GP_FCOR':detrender.data.unmasked_flux-tr_position+medsig(tr_position)[0],
                 'POS_CORR1':detrender.data.unmasked_inputs[:,1],
                 'POS_CORR2':detrender.data.unmasked_inputs[:,2]
                 })

    # lc2 = lc.copy()
    # thismask = np.copy(detrender.data.mask)

    # lc2['GP_TIME'] = np.copy(lc2[flux_type])
    # lc2['GP_TIME'][thismask] = tr_time

    # lc2['GP_POS'] = np.copy(lc2[flux_type])
    # lc2['GP_POS'][thismask] = tr_position - medsig(tr_position)[0]

    # lc2['GP_FCOR'] = np.copy(lc2[flux_type]) - lc2['GP_POS'] 
    # lc2['GP_FCOR'][lc2['GP_FCOR']==0] = np.nan

    # lc2['OUTLIER'] = np.zeros_like(lc2[flux_type]).astype('bool')
    # lc2['OUTLIER'][thismask] = detrender.data.outlier_mask.astype('bool')

    # lc2['MASK'] = np.zeros_like(lc2[flux_type]).astype('bool')
    # lc2['MASK'][thismask] = 1

    return lc2

###----------------------------------------------
###----------------------------------------------

def get_pixel_ben(ffi_name,ra,dec):
    '''Find which pixel a given star lands on at a given FFI'''

    hdre = fitsio.read_header(ffi_name,ext=ext)

    w = astropy.wcs.WCS(hdre)

    x, y = w.all_world2pix(ra,dec,0,quiet=False,adaptive=True,maxiter=100)
    
    return (x,y)

###----------------------------------------------
###----------------------------------------------

def get_mod_out(channel):
    tab = Table.read('mod_out.csv')
    index = np.where(tab['Channel']==channel)
    mod, out = tab['Mod'][index], tab['Out'][index]
    return int(mod), int(out)

###----------------------------------------------
###----------------------------------------------

def get_pixel_csv(epic,csv_file='C05_smear.csv'):
    '''Find which pixel a given star lands on at a given FFI'''

    tab = Table.read(csv_file)
    index = np.where(tab['EPIC']==int(epic))

    ra, dec, channel, col, row = float(tab['RA'][index]), float(tab['Dec'][index]), int(tab['Ch'][index]),\
        int(tab['Col'][index]), int(tab['Row'][index])
    mod, out = get_mod_out(channel)

    return ra, dec, channel, mod, out, col, row 

###----------------------------------------------
###----------------------------------------------

def get_pixel_mast(ra,dec,season):
    '''Use MAST to get a pixel position'''

    darcsec = 8.0

    cd1_1 = 0.000702794927969
    cd1_2 = -0.000853190160515
    cd2_1 = -0.000853190160515
    cd2_2 = -0.000702794927969
    cd = np.array([[cd1_1,cd1_2],[cd2_1,cd2_2]])
    cd = np.linalg.inv(cd)

    #   coordinate limits

    x1 = 1.0e30
    x2 = x1
    darcsec /= 3600.0
    ra1 = ra - darcsec / 15.0 / np.cos(dec * np.pi / 180)
    ra2 = ra + darcsec / 15.0 / np.cos(dec * np.pi / 180)
    dec1 = dec - darcsec
    dec2 = dec + darcsec

    # build mast query

    url  = 'http://archive.stsci.edu/kepler/kepler_fov/search.php?'
    url += 'action=Search'
    url += '&kic_degree_ra=' + str(ra1) + '..' + str(ra2)
    url += '&kic_dec=' + str(dec1) + '..' + str(dec2)
    url += '&max_records=100'
    url += '&verb=3'
    url += '&outputformat=JSON'

    data = json.loads(urllib.urlopen(url).read()[1:-1])

    row = data['Row_%d' % season]
    column = data['Column_%d' % season]
    channel = data['Channel_%d' % season]
    module = data['Module_%d' % season]
    output = data['Output_%d' % season]

    return int(row),int(column),int(channel),int(module),int(output)

###----------------------------------------------
###----------------------------------------------

def get_smear_file(quarter,mod,out,ddir='/kepler/kepler/smear/'):
    cads = get_num(quarter)
    fname = '%skplr%02d%1d-%d_coll.fits.gz' % (ddir,mod,out,cads)

    return fname 

###----------------------------------------------
###----------------------------------------------

def get_smear_file_k2(campaign,mod,out,ddir='/kepler/kepler2/K2/'):
    fname = '%sC%02d/collateral/ktwo%02d%1d-c%02d_coll.fits' % (ddir,int(campaign),mod,out,int(campaign))

    return fname 

###----------------------------------------------
###----------------------------------------------

def load_smear(smear_file):
    '''Return a dictionary containing the smear data and times'''
    f = fitsio.FITS(smear_file,'r')

    smear_flux = f[5]['SMEAR_FLUX'][:]
    vsmear_flux = f[3]['VSMEAR_FLUX'][:]
    time = f[5]['TIME_MJD'][:]
    cads = f[3]['CADENCENO'][:]

    smear_flux[~np.isfinite(smear_flux)] = 0 
    vsmear_flux[~np.isfinite(vsmear_flux)] = 0

    mask = np.isfinite(time)

    smear = {'smear_flux':smear_flux,
             'vsmear_flux':vsmear_flux,
             'MJD':time,
             'cads':cads}

    f.close()

    return smear

###----------------------------------------------
###----------------------------------------------

def load_output(output_file):
    '''Return a dictionary containing the smear data and times'''
    tab = Table.read(output_file)
    bjd, cad, flux = tab['BJD'], tab['CAD'], tab['FLUX']

    lc = Table({'BJD':bjd,
                'CAD': cad,
                'GP_FCOR': flux,
                'SAP_FLUX':flux
                })

    return lc

###----------------------------------------------
###----------------------------------------------


def mjd2bjd(times,ra,dec):
    '''Convert mjd to bjd - from 
    https://mail.scipy.org/pipermail/astropy/2014-April/003148.html'''
    t = Time(times, format='mjd')
    eph = jplephem.Ephemeris(de423)
    src_vec = coords.ICRS(ra=ra, dec=dec,
            unit=(units.degree, units.degree),
            distance=coords.Distance(1, units.km))
    barycenter_earthmoon = eph.position('earthmoon', t.tdb.jd)
    moonvector = eph.position('moon', t.tdb.jd)
    pos_earth = (barycenter_earthmoon - moonvector * eph.earth_share)*units.km
    corr = np.dot(pos_earth.T, src_vec.cartesian.xyz)/const.c

    dt = TimeDelta(corr, scale='tdb', format='jd')
    bjd = times+dt.value
    return bjd

###----------------------------------------------
###----------------------------------------------

def get_centroid_series(smear,x0,x1,col=None):
    '''Get the position of the star as a function of time'''

    if col is not None:
        starflux = smear[col]
    else:
        starflux = smear['smear_flux'] + smear['vsmear_flux']

    pixels = np.arange(starflux.shape[1])
    starposes = np.zeros(smear['MJD'].shape[0])
    for j in range(starflux.shape[0]):
        try:
            starposes[j] = centroid(starflux[j,x0:x1],pixels[x0:x1])
        except:
            starposes[j] = np.nan

    return starposes

###----------------------------------------------
###----------------------------------------------

def get_pixels(smear,x0,x1):
    '''Get the position of the star as a function of time'''

    starflux = np.concatenate((smear['smear_flux'][:,x0:x1],smear['vsmear_flux'][:,x0:x1]),
        axis=1)

    return Table(starflux)

###----------------------------------------------
###----------------------------------------------

def get_background(smear,col=None,cutoff=25):
    '''Get the background flux as a function of time, by smoothing a mean of
    low-flux columns'''

    if col is not None:
        starflux = smear[col]
    else:
        starflux = smear['smear_flux'] + smear['vsmear_flux']

    mean_flux = np.mean(starflux,axis=0)
    back_cols = (mean_flux<np.percentile(mean_flux,cutoff))

    raw_background = np.mean(starflux[:,back_cols],axis=1)
    raw_background[raw_background<0.05*np.median(raw_background)] = np.nan

    background = savgol_filter(raw_background,7,5)#NIF(raw_background,250,11)

    return background

###----------------------------------------------
###----------------------------------------------

def gpfilt_1(flux,time,scale,yerr=0.01):
    '''Do a 1D squared exponential GP fit'''

    keep = np.isfinite(flux)

    # initialize a kernel
    k1 = 1.0*kernels.ExpSquaredKernel(1.0)
    gp = gp = george.GP(k1)
    gp.compute(time[keep],yerr=yerr)

    # do a first prediction
    ypred, ycov = gp.predict(flux[keep],time)
    err_pred = np.sqrt(np.diag(ycov)+yerr**2)

    return ypred, err_pred, gp

###----------------------------------------------
###----------------------------------------------

def gp_negll(par, X, y, gp):
    gp.kernel.pars[0] = 10**par[0]
    gp.compute(X, yerr = 10**par[1], sort=False)
    return -gp.lnlikelihood(y)

###----------------------------------------------
###----------------------------------------------

def GP_train_bound(x, y, gp, par_in, bounds):
    ''' 
    Max likelihood optimization of GP hyper-parameters within bounds.
    '''
    args = (x, y, gp)
    res = op.fmin_l_bfgs_b(gp_negll, par_in, \
                            args = args, approx_grad = True, \
                            bounds = bounds, disp = 0)
    res[0][0] = 10**res[0][0]
    return res[0]

###----------------------------------------------
###----------------------------------------------

def gpfilt_iter(flux,time,scale,yerr=0.1,iters=5,nsig=3.0,bounds=((-1,2),(-3,None))):
    '''Do a 1D squared exponential GP fit'''
    mask = np.isfinite(flux)

    work = np.copy(flux[mask])
    t = np.copy(time[mask])

    for j in range(iters):
        work, err_pred, gp = gpfilt_1(work,t,scale,yerr=yerr)
        out = np.isnan(work) + myabs((flux-work)>(nsig*err_pred))
        work[out] = np.nan
        scale, yerr = 10**(GP_train_bound(t[~out],work[~out],gp,
            (np.log10(scale),np.log10(yerr)),bounds))
    print scale, yerr
    work, err_pred, gp = gpfilt_1(work,t,scale,yerr=yerr)
    
    return work

###----------------------------------------------
###----------------------------------------------

def extract_lc(smear,starposes,background,width,col=None):
    '''Extract a light curve, given positions, backgrounds and aperture'''

    if col is not None:
        starflux = smear[col]
    else:
        starflux = smear['smear_flux'] + smear['vsmear_flux']

    poses = np.arange(starflux.shape[1])

    flux = np.zeros(smear['MJD'].shape[0])

    for j in range(starflux.shape[0]):
        if ~np.isfinite(starposes[j]):
            # check it's a good frame
            flux[j] = np.nan
            continue

        # make a psf
        # psf = cosbell(poses,np.round(starposes[j]*4)/4.,width,taper,subsample=8)
        psf = supergaussian_bin(poses,starposes[j],width)
        psf /= np.nanmax(psf)

        #extract a flux
        flux[j] = np.sum((psf)*starflux[j,:]) - np.sum(psf)*background[j]

    return flux

###----------------------------------------------
###----------------------------------------------


def do_target(name,quarter,cat_file='kepler_inputs.csv',out_dir = 'kepler_smear/',
    cbvdir = '/kepler/kepler/FITS/',smear_type=None,do_plot=False,
    gap_file = '/users/popeb/K2LC/ben/kepler_bad_cads.csv',sig_clip=True):
    
    tab = Table.read(cat_file)
    smear_name = lambda s: '' if s is None else str(s)
    index = np.where(tab['Name']==name)
    rah, dech = tab['RA'][index][:].data[0], tab['Dec'][index][:].data[0]
    ra, dec = sex2dec(rah,dech)
    print '\nExtracting smear light curve for %s , quarter %d' % (name,quarter) 
    print 'RA %f, Dec %f' % (ra, dec)

    # get star location 
    season = np.mod(quarter+2,4)

    print '\nQuerying MAST for star position, season %d' % season

    row,col,channel,mod,out = get_pixel_mast(ra,dec,season)
    
    print 'Found star at mod.out %d.%d, channel %d, row %d, column %d' % \
        (mod,out,channel,row,col)

    col -= 12

    # load the smear data 

    print '\nLoading smear data'
    
    fname = get_smear_file(quarter,mod,out)

    smear = load_smear(fname)
    print 'Smear data loaded'

    ## now do MJD - bjd correction
    mjd = np.copy(smear['MJD'][np.isfinite(smear['MJD'])])
    bjd = np.copy(smear['MJD'])
    bjd[np.isfinite(smear['MJD'])] = mjd2bjd(mjd,ra,dec)
    t0 = np.nanmin(bjd)

    print 'Corrected MJD to BJD'

    if do_plot:
        plt.clf()
        plt.plot(smear['smear_flux'].mean(axis=0),'b')
        plt.plot(smear['vsmear_flux'].mean(axis=0),'g')
        plt.axvline(col,color='r')
        plt.xlabel('pix')
        plt.ylabel('Flux (counts)')
        plt.title('Smear Profile, mod.out %d.%d, Q%d' % (mod,out,quarter))
        plt.savefig('%sprofile_%d%d_q%d.png' % (out_dir,mod,out,quarter))
        print 'Saved smear profile to %sprofile_%d%d_q%d.png' % (out_dir,mod,out,quarter)

    # find the star 
    print '\nFitting star centroids'

    starposes = get_centroid_series(smear,col-5,col+5,col=smear_type)

    if do_plot:
        plt.clf()
        plt.plot(bjd-t0,starposes-np.nanmedian(starposes),'.k')
        plt.xlabel('BJD - %f' % t0)
        plt.ylabel(r'$\Delta$ pix')
        plt.title('%s position, Q%d' % (name,quarter))
        plt.savefig('%sstarpos_%s_q%d.png' % (out_dir,name,quarter))
        print 'Saved star positions to %sstarpos_%s_q%d.png' % (out_dir,name,quarter)

    # get background flux 
    print '\nExtracting background flux'

    background = get_background(smear,col=smear_type)

    if do_plot:
        plt.clf()
        plt.plot(bjd-t0,background)
        plt.xlabel('BJD - %f' % t0)
        plt.ylabel('Flux (counts/pix)')
        plt.title('Background flux, mod.out %d.%d, Q%d' % (mod,out,quarter))
        plt.savefig('%sbackground_%d%d_q%d.png' % (out_dir,mod,out,quarter))
        print 'Saved background to %sbackground_%d%d_q%d.png' % (out_dir,mod,out,quarter)

    # now extract a lightcurve for several different apertures
    print '\nExtracting aperture photometry'

    widths = [1.5,2,3,4,5]
    # tapers = [2,3,3,3,4]

    flux = np.zeros((5,bjd.shape[0]))
    sigs = np.zeros(5)

    for j in range(5):
        flux[j,:] = extract_lc(smear,starposes,background,
            widths[j],col=smear_type)
        sigs[j] = cdpp(bjd,flux[j,:])

    best_aperture = np.argmin(sigs)

    print 'Best aperture:',best_aperture

    if do_plot:
        plt.clf()
        for j in range(5):
            plt.plot(bjd-t0,flux[j,:],'.')
        plt.xlabel('BJD -%f' % t0)
        plt.ylabel('Flux (counts)')
        plt.title('%s Light curves, Q%d' % (name,quarter))
        plt.savefig('%slcs_%s_q%d.png' % (out_dir,name,quarter))
        print 'Saved light curves to %slcs_%s_q%d.png' % (out_dir,name,quarter)

        plt.clf()
        plt.plot(bjd-t0,flux[best_aperture,:]/
            (medsig(flux[best_aperture,:])[0]),'.k')
        plt.xlabel('BJD -%f' % t0)
        plt.ylabel('Relative Flux')
        plt.title('%s Optimal Light curve, Q%d' % (name,quarter))
        plt.savefig('%slc_best_%s_q%d.png' % (out_dir,name,quarter))
        print 'Saved best light curve to %slc_best_%s_q%d.png' % (out_dir,name,quarter)


    # get quality flags
    try:
        nearby_tab = my_MASTRADec(ra,dec,tel='kepler',quarter=quarter)

        quality = nearby_tab['SAP_QUALITY']
    except:
        print 'Expanding search!'
        nearby_tab = my_MASTRADec(ra,dec,tel='kepler',quarter=quarter,darcsec=50000)

        quality = nearby_tab['SAP_QUALITY']

    lc = Table({'BJD':bjd,
                'MJD':smear['MJD'],
                'CAD': smear['cads'],
                'FLUX1': flux[0,:],
                'FLUX2': flux[1,:],
                'FLUX3': flux[2,:],
                'FLUX4': flux[3,:],
                'FLUX5': flux[4,:],
                'FLUX' : flux[best_aperture,:],
                'BACKGROUND':background,
                'STARPOS':starposes,
                'QUALITY':quality
                })

    # censor bad cadences
    for bit in [20,21]:
        bad = (quality & (2**bit)) == (2**bit)

        for key in lc.keys():
            if 'FLUX' in key:
                lc[key][bad] = np.nan         

    censoredlc = censor_bad_cads(lc,quarter,gap_file=gap_file)

    if do_plot:
        plt.clf()
        plt.plot(lc['BJD'],lc['FLUX'],'.r')
        plt.plot(censoredlc['BJD'],censoredlc['FLUX'],'.k')
        plt.xlabel('BJD')
        plt.ylabel('FLUX')
        plt.title('')
        plt.savefig('%slc_bads_%s_q%d.png' % (out_dir,name,quarter))
        print 'Saved corrected light curve to %slc_bads_%s_q%d.png' % (out_dir,name,quarter)

    lc = censoredlc.copy()
    # now do a jump correction

    # jumplc = smear_jump(lc,name,do_plot=True)
    # if do_plot:
    #   plt.savefig('%slc_jumps_%s_q%d.png' % (out_dir,name,quarter))
    #   print 'Saved corrected light curve to %slc_jumps_%s_q%d.png' % (out_dir,name,quarter)
    #   plt.clf()

    # cbv detrending

    cbvfile = glob('%s*q%02d*.fits' % (cbvdir,quarter))[0] # this should be unique
    
    # do this for both 4 and 8 cbvs just to check diversity
    flux4s = np.zeros((5,bjd.shape[0]))
    flux8s = np.zeros((5,bjd.shape[0]))

    for aperture in range(5):
        print 'Doing aperture', aperture
        dummy = lc.copy()
        dummy['FLUX'] = dummy['FLUX%d' % (aperture+1)]

        ### First do 4 CBVs

        cbtime, cbcadence, cbraw_flux, flux_cbv4, cbweights = cbv.correct_smear(dummy, 
            cbvfile, name, quarter, mod,out, nB = 4, outfile = None, 
            exclude_func = None, exclude_func_par = None, doPlot = True)

        if sig_clip:
            dummy = lc.copy()
            smooth4 = NIF(flux_cbv4,101,15)
            white4 = dummy['FLUX'] - smooth4
            mm, ss = medsig(white4)
            outliers = (np.abs(white4)> (4*ss))
            if np.sum(outliers) < (0.2*ncad):
                dummy['FLUX'][outliers] = np.nan
                cbtime, cbcadence, cbraw_flux, flux_cbv4, cbweights = cbv.correct_smear(dummy, 
                    cbvfile, name, quarter, mod,out, nB = 4, outfile = None, 
                    exclude_func = None, exclude_func_par = None, doPlot = True)

        if do_plot:
            plt.savefig('%slc_corr4_%s_q%d.png' % (out_dir,name,quarter))
            print 'Saved corrected light curve to %slc_corr4_%s_q%d.png' % (out_dir,name,quarter)
            plt.clf()

        ### Now do 8 CBVs

        dummy = lc.copy()
        dummy['FLUX'] = dummy['FLUX%d' % (aperture +1)]
        cbtime, cbcadence, cbraw_flux, flux_cbv8, cbweights = cbv.correct_smear(dummy, 
            cbvfile, name, quarter, mod,out, nB = 8, outfile = None, 
            exclude_func = None, exclude_func_par = None, doPlot = True)

        if sig_clip:
            dummy = lc.copy()
            smooth8 = NIF(flux_cbv8,101,15)
            white8 = dummy['FLUX'] - smooth8
            mm, ss = medsig(white8)
            outliers = (np.abs(white8)> (4*ss))
            if np.sum(outliers) < (0.2*ncad):
                dummy['FLUX'][outliers] = np.nan
                cbtime, cbcadence, cbraw_flux, flux_cbv8, cbweights = cbv.correct_smear(dummy, 
                    cbvfile, name, quarter, mod,out, nB = 8, outfile = None, 
                    exclude_func = None, exclude_func_par = None, doPlot = True)

        if do_plot:
            plt.savefig('%slc_corr8_%s_q%d.png' % (out_dir,name,quarter))
            print 'Saved corrected light curve to %slc_corr8_%s_q%d.png' % (out_dir,name,quarter)
            plt.clf()

        # now save the lightcurve
        flux4s[aperture,:] = flux_cbv4
        flux8s[aperture,:] = flux_cbv8

        lc['FLUX%d_CORR_4' % aperture] = flux_cbv4
        lc['FLUX%d_CORR_8' % aperture] = flux_cbv8

    sigs4 = np.zeros(5)
    sigs8 = np.zeros(5)

    for j in range(5):
        # mm, ss = medsig(flux4s[j,:])
        sigs4[j] = cdpp(lc['BJD'],flux4s[j,:])#ss/mm

    best_aperture4 = np.argmin(sigs4)

    for j in range(5):
        # mm, ss = medsig(flux8s[j,:])
        sigs8[j] = cdpp(lc['BJD'],flux8s[j,:])

    best_aperture8 = np.argmin(sigs8)

    print 'Best Aperture (4 CBVs):', best_aperture4
    print 'Best Aperture (8 CBVs):', best_aperture8

    lc['FLUX_CORR_4'] = flux4s[best_aperture4]
    lc['FLUX_CORR_8'] = flux8s[best_aperture8]

    return lc 

###----------------------------------------------
###----------------------------------------------

def do_target_k2(name,campaign,cat_file='k2_inputs.csv',out_dir = 'k2_smear/',
    smear_type=None,do_plot=False,do_pixels=False):
    
    tab = Table.read(cat_file)
    smear_name = lambda s: '' if s is None else str(s)
    index = np.where(tab['EPIC']==int(name))
    rah, dech = tab['RA'][index][:].data[0], tab['Dec'][index][:].data[0]
    try:
        ra, dec = sex2dec(rah,dech)
    except:
        ra, dec = rah, dech
    print '\nExtracting smear light curve for %s , Campaign %s' % (name,campaign) 

    # get star location 

    ra, dec, channel, mod, out, col, row = get_pixel_csv(name,csv_file=cat_file)
    print 'RA %f, Dec %f' % (ra, dec)
    
    print 'Found star at mod.out %d.%d, channel %d, row %d, column %d' % \
        (mod,out,channel,row,col)

    col -= 12

    # load the smear data 

    print '\nLoading smear data'
    
    fname = get_smear_file_k2(int(campaign[-1]),mod,out)

    smear = load_smear(fname)
    print 'Smear data loaded'

    ## now do MJD - bjd correction
    mjd = np.copy(smear['MJD'][np.isfinite(smear['MJD'])])
    bjd = np.copy(smear['MJD'])
    bjd[np.isfinite(smear['MJD'])] = mjd2bjd(mjd,ra,dec)
    t0 = np.nanmin(bjd)

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

    if do_pixels:
        pixels = get_pixels(smear,np.nanmedian(starposes)-5,np.nanmedian(starposes)+5)
        pixels.write('%s%s_%s_pixels_%d%d.fits' % (out_dir,campaign,name,mod,out),overwrite=True)
        print '%d Pixels saved to %s%s_%s_pixels_%d%d.fits' % (len(pixels.keys()),out_dir,campaign,name,mod,out)

    if do_plot:
        plt.clf()
        plt.plot(bjd-t0,starposes-np.nanmedian(starposes),'.k')
        plt.xlabel('BJD - %f' % t0)
        plt.ylabel(r'$\Delta$ pix')
        plt.title('%s position' % (name))
        plt.savefig('%s%s_starpos.png' % (out_dir,name))
        print 'Saved star positions to %s%s_starpos.png' % (out_dir,name)

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

    widths = [1.5,2,3,4,5]

    flux = np.zeros((5,bjd.shape[0]))
    sigs = np.zeros(5)

    for j in range(5):
        flux[j,:] = extract_lc(smear,starposes,background,
            widths[j],col=smear_type)
        sigs[j] = cdpp(bjd,flux[j,:])
        if ~np.isfinite(sigs[j]):  
            sigs[j] = 1e10

    best_aperture = np.argmin(sigs)

    print 'Best aperture:',best_aperture

    if do_plot:
        plt.clf()
        for j in range(5):
            plt.plot(bjd-t0,flux[j,:],'.')
        plt.xlabel('BJD -%f' % t0)
        plt.ylabel('Flux (counts)')
        plt.title('%s Light curves' % (name))
        plt.savefig('%s%s_lcs.png' % (out_dir,name))
        print 'Saved light curves to %s%s_lcs.png' % (out_dir,name)

        plt.clf()
        plt.plot(bjd-t0,flux[best_aperture,:]/
            (medsig(flux[best_aperture,:])[0]),'.k')
        plt.xlabel('BJD -%f' % t0)
        plt.ylabel('Relative Flux')
        plt.title('%s Optimal Light curve' % (name))
        plt.savefig('%s%s_lc_best.png' % (out_dir,name))
        print 'Saved best light curve to %s%s_lc_best.png' % (out_dir,name)
            # query mast for the nearest object

    try:
        nearby_tab = my_MASTRADec(ra,dec)

        x, y, quality = nearby_tab['POS_CORR1'], nearby_tab['POS_CORR2'],\
            nearby_tab['SAP_QUALITY']

        print '\nSet metadata equal to nearby star'
    except:
        print 'Failed to download a nearby light curve'
        try:
            print 'Expanding search!'
            nearby_tab = my_MASTRADec(ra,dec,darcsec=50000)
            x, y, quality = nearby_tab['POS_CORR1'], nearby_tab['POS_CORR2'],\
                nearby_tab['SAP_QUALITY']
            print '\nSet metadata equal to nearby star'
        except:
            print 'Failed totally to find a nearby light curve'
            quality = ~np.isfinite(flux[0,:])
            x, y = np.zeros_like(bjd), np.zeros_like(bjd)


    lc = Table({'BJD':bjd,
                'TIME':bjd,
                'MJD':smear['MJD'],
                'CADENCENO': smear['cads'],
                'FLUX1': flux[0,:],
                'FLUX2': flux[1,:],
                'FLUX3': flux[2,:],
                'FLUX4': flux[3,:],
                'FLUX5': flux[4,:],
                'SAP_FLUX' : flux[best_aperture,:],
                'SAP_FLUX_ERR' : flux[np.argmin(sigs),:],
                'SAP_QUALITY':quality,
                'POS_CORR1': x,
                'POS_CORR2': y,
                'BACKGROUND':background,
                'STARPOS':starposes,
                })

    return lc 

###----------------------------------------------
###----------------------------------------------

smear_name = lambda s: '' if s is None else str(s)

###----------------------------------------------
###----------------------------------------------

def stitch_lcs(out_dir,name,smear_type,do_plot=True):
    flux = np.array([])
    flux4 = np.array([])
    flux8 = np.array([])
    allf4s = np.zeros((5,0))
    allf8s = np.zeros((5,0))
    time = np.array([])
    smooth = np.array([])
    quarterlist = np.array([])
    medflux = np.array([])
    cadence = np.array([])
    background = np.array([])

    for j in range(18):
        try:
            lcfile = '%s%s_smear_q%d%s.csv' % (out_dir,name,j,smear_name(smear_type))
            lc = Table.read(lcfile)

            thisflux = lc['FLUX']
            thisflux4 = lc['FLUX_CORR_4']
            thisflux8 = lc['FLUX_CORR_8']
            thesecads = lc['CAD']
            ncad = np.size(thesecads)
            f4s  = np.zeros((5,ncad))
            f8s = np.zeros((5,ncad))

            for k in range(5): # do all the corrected fluxes
                f4s[k] = lc['FLUX%d_CORR_4' % k]
                f8s[k] = lc['FLUX%d_CORR_8' % k]
                f4s[k] /= medsig(f4s[k])[0]
                f8s[k] /= medsig(f8s[k])[0]

            medflux = np.append(medflux,np.ones_like(thisflux)*medsig(thisflux)[0])
            thisflux /= medsig(thisflux)[0]
            thisflux4 /= medsig(thisflux4)[0]
            thisflux8 /= medsig(thisflux8)[0]
            thistime = lc['BJD']
            thissmooth = NIF(thisflux8,101,11)
            thisq = j*np.ones_like(thisflux)

            flux = np.concatenate((flux,thisflux))
            flux4 = np.concatenate((flux4,thisflux4))
            flux8 = np.concatenate((flux8,thisflux8))
            allf4s = np.concatenate((allf4s,f4s),axis=1)
            allf8s = np.concatenate((allf8s,f8s),axis=1)
            time = np.concatenate((time,thistime))
            smooth = np.concatenate((smooth,thissmooth))
            quarterlist = np.concatenate((quarterlist,thisq))
            background = np.concatenate((background,lc['BACKGROUND']))
            cadence = np.concatenate((cadence,thesecads))
        except:
            print 'Failed on Quarter',j
            continue

    # flux *= medsig(medflux)[0]
    # flux4 *= medsig(medflux)[0]
    # flux8 *= medsig(medflux)[0]
    # smooth *= medsig(medflux)[0]

    newlc = Table({'BJD':time,
                   'CADENCENO':cadence,
                   'FLUX':flux,
                   'FLUX_CORR_4':flux4,
                   'FLUX_CORR_8':flux8,
                   'FILT':smooth,
                   'QUARTER':quarterlist,
                   'BACKGROUND':background,
                   'MEDFLUX':medflux
                })

    for j in range(5):
        newlc['FLUX%d_CORR_4' % j] = allf4s[j,:]
        newlc['FLUX%d_CORR_8' % j] = allf8s[j,:]

    if do_plot:
        plt.clf()
        plt.plot(time,flux8,'.k')
        plt.xlabel('BJD')
        plt.ylabel('Flux')
        plt.title('%s Optimal Light curve' % (name))
        plt.savefig('%slc_best_%s.png' % (out_dir,name))
        print 'Saved best light curve to %slc_best_%s.png' % (out_dir,name)
        
    return newlc

def stitch_combine(out_dir,name,do_plot=True,thresh=2.5):
    ''' you can tell the provenance by the new column SMEAR_TYPE which is 2 for both,
    0 for masked, 1 for virtual'''

    lc_masked = Table.read('%s%s_smear_full%s.csv' % (out_dir,name,smear_name('smear_flux')))
    lc_virtual = Table.read('%s%s_smear_full%s.csv' % (out_dir,name,smear_name('vsmear_flux')))
    lc_tot = Table.read('%s%s_smear_full.csv' % (out_dir,name))

    lc_masked['SMEAR_TYPE'] = np.zeros_like(lc_masked['FLUX'])
    lc_virtual['SMEAR_TYPE'] = np.ones_like(lc_virtual['FLUX'])

    mask_sigs = []
    virtual_sigs = []

    for q in range(17):
        m = lc_tot['QUARTER'] == q

        mmed, msig = medsig(lc_masked['FLUX_CORR_8'][m])
        vmed, vsig = medsig(lc_virtual['FLUX_CORR_8'][m])

        mask_sigs.append(msig/mmed)
        virtual_sigs.append(vsig/vmed)

    mask_sigs = np.array(mask_sigs)
    virtual_sigs = np.array(virtual_sigs)

    medmasksigs = medsig(mask_sigs)[0]
    medvirtsigs = medsig(virtual_sigs)[0]

    # now go back over and pick good ones

    dummy = lc_tot.copy()
    dummy['SMEAR_TYPE'] = 2.*np.ones_like(dummy['FLUX'])

    for q in range(17):
        try:
            m = lc_tot['QUARTER'] == q
            if mask_sigs[q] >= thresh*(virtual_sigs[q]):
                dummy[m] = lc_virtual[m]
                print 'Using masked smear for quarter',q
            elif virtual_sigs[q] >= thresh*(mask_sigs[q]):
                dummy[m] = lc_masked[m]
                print 'Using virtual smear for quarter',q
        except:
            pass

    return dummy 