import pylab, numpy, pyfits
import os
from pylab import *
from matplotlib import *
from numpy import *
import sys, urllib, time, re, json
import fitsio
import wget
from PIL import Image

from astroquery.simbad import Simbad 
from astroquery.skyview import SkyView
from astropy.table import Table

skygroup = ''; channel = ''; module = ''; output = ''; column = ''; row = ''

def simbad_search(ra,dec):
    return Simbad.query_region(ra+' '+dec)

def get_postage_stamp_jpeg(epic,ra,dec,surveys=['DSS2 Red','DSS2 Blue', 'DSS2 IR']):
    path = SkyView.get_images(position=ra+' '+dec,survey=surveys)

    sz = path[0][0].data.shape[0]
    r, g, b = path[0][0].data, path[1][0].data, path[2][0].data
    r, g, b = r/r.max(), g/g.max(), b/b.max()

    rgb_uint8 = (np.dstack((r,g,b)) * 255.999) .astype(np.uint8)

    img = Image.fromarray(rgb_uint8)

    img.save('EPIC%d_DSS_im.jpg' % epic)

    print 'Saved image to EPIC%d_DSS_im.jpg' % epic


def MAST_EPIC(id):
    '''Hacked from Tom Barclay's kepffi.py'''

    global skygroup, column, row
    season = 0
    id = str(id)
    
# build mast query

    url  = 'http://archive.stsci.edu/k2/data_search/search.php?'
    url += 'action=Search'
    url += '&ktc_k2_id=' + id
    url += '&max_records=100'
    url += '&verb=3'
    url += '&outputformat=JSON'

# retrieve results from MAST

    out = ''
    try:
        lines = urllib.urlopen(url)
        raw = lines.read()
        data = json.loads(raw)[0]
        # for line in lines:
        #     line = line.strip()
        #     if (len(line) > 0 and 
        #         'Kepler' not in line and 
        #         'integer' not in line and
        #         'no rows found' not in line):
        #         out = line.split(',')
        # if len(out) > 0:
        #     for jj in [4,13,15,19,21,25,27,29,33,35,37]:
        #         try:
        #             out[jj] = float(out[jj])
        #         except:
        #             pass
        #     data = {'EPIC':out[0],
        #             'RA': out[1],
        #             'Dec':out[2],
        #             'Avail':out[3],
        #             'Kp_mag':(out[4]),
        #             'HIP':out[5],
        #             'TYC':out[6],
        #             'UCAC':out[7],
        #             '2MASS':out[8],
        #             'pmra':(out[13]),
        #             'pmdec':(out[15]),
        #             'Bmag':(out[19]),
        #             'Vmag':(out[21]),
        #             'gmag':(out[25]),
        #             'rmag':(out[27]),
        #             'imag':(out[29]),
        #             'Jmag':(out[33]),
        #             'Hmag':(out[35]),
        #             'Kmag':(out[37])
        #             }

    except:
        txt = 'ERROR -- no target found with KepID %s' % id
        sys.exit(txt)

    return data#kepid,ra,dec,kepmag#,skygroup,channel,module,output,row,column

def my_MASTRADec(ra,dec,darcsec=1200.,tel='k2',quarter=None):

    global skygroup, column, row

# WCS data

    # cd1_1 = 0.000702794927969
    # cd1_2 = -0.000853190160515
    # cd2_1 = -0.000853190160515
    # cd2_2 = -0.000702794927969
    # cd = array([[cd1_1,cd1_2],[cd2_1,cd2_2]])
    # cd = linalg.inv(cd)

# coordinate limits

    x1 = 1.0e30
    x2 = x1
    darcsec /= 3600.0
    ra1 = ra - darcsec / 15.0 / cos(dec * pi / 180)
    ra2 = ra + darcsec / 15.0 / cos(dec * pi / 180)
    dec1 = dec - darcsec
    dec2 = dec + darcsec

# build mast query

    url  = 'http://archive.stsci.edu/%s/data_search/search.php?' % tel
    url += 'action=Search'
    if tel =='kepler':
        url+= '&sci_data_quarter='+str(quarter)
        url+='&ktc_target_type=LC'
    url += '&ra=' + str(ra1) + '..' + str(ra2)
    url += '&dec=' + str(dec1) + '..' + str(dec2)
    url += '&max_records=100'
    url += '&verb=3'
    # url += '&outputformat=CSV'
    url += '&outputformat=WGET_file'

# retrieve results from MAST: nearest EPIC source to supplied coordinates
    out = ''
    lines = urllib.urlopen(url)
    for line in lines:
        line = line.strip()
        if (len(line) > 0 and 'not' not in line and 
            'Kepler' not in line and 
            'integer' not in line and
            'no rows found' not in line):
            out = line.split(',')

    url = out[0][8:]

    # download only the first URL

    fname = wget.download(url)

    tab = Table.read(fname)
    try:
        os.remove(fname)
    except:
        pass

    return tab