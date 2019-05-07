#% matplotlib inline
#%config InlineBackend.figure_format = 'retina'
#%load_ext line_profiler
#%load_ext autoreload
#%autoreload 2

from __future__ import division

import numpy as np
import glob, os, json, pickle, sys
import matplotlib.pyplot as plt
import scipy.linalg as sl

import libstempo, libstempo.toasim as LT
import enterprise
from enterprise.pulsar import Pulsar
import enterprise.signals.parameter as parameter
from enterprise.signals import utils
from enterprise.signals import signal_base
from enterprise.signals import selections
from enterprise.signals.selections import Selection
from enterprise.signals import white_signals
from enterprise.signals import gp_signals
from enterprise.signals import deterministic_signals

import corner
from PTMCMCSampler.PTMCMCSampler import PTSampler as ptmcmc


import models
from models import models, model_utils

#Get the name of directory corresponding to amplitude of gwb:
loc = np.int(sys.argv[2])
A_gwb_1 = np.linspace(1e-16, 1e-15, 20)
A_gwb_2 = np.linspace(1.1e-15, 5e-15, 10)

A_gwb = np.append(A_gwb_1, A_gwb_2)

#injection_dir =	'injections/' +	'injecting_' + str(A_gwb[loc]) + '_gwb/'

#print injection_dir

realiz = np.int(sys.argv[1])

datadir = '/gpfs/home/nspol/stochastic_11yr_analysis/data/'
realiz_dir = '/gpfs/scratch/nspol/injected_timfiles/realization_' + str(realiz) + '/'
injection_dir = 'injecting_' + str(A_gwb[loc]) + '_gwb/'

def get_noise_from_enterprise(noisefile):
    
    with open(noisefile) as f:
        params = json.load(f)
    
    return params

parfiles = sorted(glob.glob(datadir + 'partim/*.par'))
timfiles = sorted(glob.glob(realiz_dir + injection_dir + '*.tim'))
noisefiles = sorted(glob.glob(datadir + 'noisefiles/*.json'))

psrs = []
for p, t in zip(parfiles, timfiles):
    psr = Pulsar(p, t, ephem='DE436')
    psrs.append(psr)
    
params = {}
for nfile in noisefiles:
    params.update(get_noise_from_enterprise(nfile))
       
pta = models.model_3a(psrs, psd='powerlaw', noisedict=params, gamma_common=13./3,
                      upper_limit=False, bayesephem=True)

outdir = '/scratch/nspol/11_yr_injection_3a/' + 'realization_' + str(realiz) + '/' + 'injection_' + str(loc) + '/'
sampler = model_utils.setup_sampler(pta, resume=True, outdir=outdir)

N = int(5e6) # one mega-sample!
x0 = np.hstack(p.sample() for p in pta.params)
sampler.sample(x0, N, AMweight=25, SCAMweight=40, DEweight=55)
