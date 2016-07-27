# keplersmear

[![Licence](http://img.shields.io/badge/license-GPLv3-blue.svg?style=flat)](http://www.gnu.org/licenses/gpl-3.0.html)
[![arXiv](http://img.shields.io/badge/arXiv-1603.09167-blue.svg?style=flat)](http://arxiv.org/abs/1510.00008)

Python package for generating light curves from Kepler and K2 collateral data. 

Using smear data, you can recover photometry of bright stars (Kp <~ 8) that were saturated or otherwise not conventionally observed by Kepler/K2. 

This is an evolving body of code with complicated dependencies. We are endeavouring to simplify this and produce a more straightforward standalone package. 

## Installation

    git clone https://github.com/benjaminpope/keplersmear/
    cd keplersmear
    python setup.py install --user

## Requires

 - NumPy, SciPy, astropy, jplephem, de423, George, SuzPyUtils, fitsio, PyKE

## Detrending requires

 - k2sc (K2), KeplerSys (nominal Kepler)

Citing
------

If you use K2SC in your research, please cite

	Pope, B. J. S.; White, T. R.; Huber, D.; Murphy, S. J.; Bedding, T. R.; Caldwell, D. A.; Sarai, A.; Aigrain, S.; Barclay, T. (MNRAS, 2016), arXiv:1510.00008

or use this ready-made BibTeX entry

	@ARTICLE{2016MNRAS.455L..36P,
	   author = {{Pope}, B.~J.~S. and {White}, T.~R. and {Huber}, D. and {Murphy}, S.~J. and 
		{Bedding}, T.~R. and {Caldwell}, D.~A. and {Sarai}, A. and {Aigrain}, S. and 
		{Barclay}, T.},
	    title = "{Photometry of very bright stars with Kepler and K2 smear data}",
	  journal = {\mnras},
	archivePrefix = "arXiv",
	   eprint = {1510.00008},
	 primaryClass = "astro-ph.IM",
	 keywords = {asteroseismology, techniques: photometric, stars: individual: HR 8500, 70 Aqr, HD 178875, stars: variables: general},
	     year = 2016,
	    month = jan,
	   volume = 455,
	    pages = {L36-L40},
	      doi = {10.1093/mnrasl/slv143},
	   adsurl = {http://adsabs.harvard.edu/abs/2016MNRAS.455L..36P},
	  adsnote = {Provided by the SAO/NASA Astrophysics Data System}
	}
