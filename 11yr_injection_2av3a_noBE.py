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

def get_noise_from_enterprise(noisefile):
    
    with open(noisefile) as f:
        params = json.load(f)
    
    return params

realiz = np.int(sys.argv[1])
loc = np.int(sys.argv[2])

A_gwb_1 = np.linspace(1e-16, 1e-15, 20)
A_gwb_2 = np.linspace(1.1e-15, 5e-15, 10)

A_gwb = np.append(A_gwb_1, A_gwb_2)

datadir = '/gpfs/home/nspol/stochastic_11yr_analysis/data/'
realiz_dir = '/gpfs/scratch/nspol/real_injected_timfiles/realization_' + str(realiz) + '/'
injection_dir = 'injecting_' + str(A_gwb[loc]) + '_gwb/'

parfiles = sorted(glob.glob(datadir + 'partim/*.par'))
timfiles = sorted(glob.glob(realiz_dir + injection_dir + '*.tim'))
#noisefiles = sorted(glob.glob(datadir + 'noisefiles/*.json'))
psrlist = np.loadtxt("psrlist.txt", dtype = 'string')

noisefiles_path = datadir + 'noisefiles/'

psrs = []
params = {}

for p, t in zip(parfiles, timfiles):
    psr = Pulsar(p, t, ephem='DE436')
    #Check if this pulsar has baseline > 3 years
    psr_name = psr.name
    if psr_name in psrlist:
        
        psrs.append(psr)
        full_path = noisefiles_path + '*' + str(psr_name) + '*'
        noisefile = glob.glob(full_path)
        
        params.update(get_noise_from_enterprise(noisefile[0]))
        
    else:
        continue

#Initialize dictionary with two keys for model 2a and 3a
pta = dict.fromkeys(np.arange(0, 2))

pta[0] = models.model_2a(psrs, psd = 'powerlaw', noisedict = params, gamma_common = 13./3,
                      upper_limit = False, bayesephem = False)

pta[1] = models.model_3a(psrs, psd = 'powerlaw', noisedict = params, gamma_common = 13./3, 
                         upper_limit = False, bayesephem = False)

super_model = model_utils.HyperModel(pta)

outdir = '/scratch/nspol/real_injected_results/2av3a_wo_BE/' + 'realization_' + str(realiz) + '/' + 'injection_' + str(loc)
sampler = super_model.setup_sampler(resume=True, outdir=outdir)

N = int(5e6) # one mega-sample!
#x0 = np.hstack(p.sample() for p in pta.params)
x0 = super_model.initial_sample()
sampler.sample(x0, N, AMweight=25, SCAMweight=40, DEweight=55)