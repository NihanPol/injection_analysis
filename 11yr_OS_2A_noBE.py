from __future__ import division

import pandas as pd, numpy as np, glob
import matplotlib.pyplot as pl

import numpy as np, pandas as pd
import glob, os as OS, json, sys
import matplotlib.pyplot as plt
import scipy.linalg as sl
import pandas as pd

import acor

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

#Import optimal statistic stuff
from enterprise_extensions import model_utils
from enterprise_extensions import models
from enterprise_extensions.frequentist import optimal_statistic

import corner
from PTMCMCSampler.PTMCMCSampler import PTSampler as ptmcmc

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

chain_dir = '/scratch/nspol/real_injected_results/2a_wo_BE/' + \
'realization_' + str(realiz) + '/' + 'injection_' + str(loc) + '/'

chain = np.loadtxt(chain_dir + '/chain_1.txt')
burn = int(0.25*chain.shape[0])
pars = np.loadtxt(chain_dir + '/pars.txt', dtype=np.unicode_)

os = optimal_statistic.OptimalStatistic(psrs, gamma_common = 13./3.,  orf = 'hd', bayesephem = False, noisedict = params)
#os.pta.set_default_params(params)

opt, snr = os.compute_noise_marginalized_os(chain[burn:, :], N = int(1e4))

outdir = '/gpfs/scratch/nspol/os_output/real_data/model_2a_noBE/'
dirname = "realization_" + str(realiz) + '/'

op_fname = outdir + dirname + 'amp_' + str(loc) + '.npz'

if not OS.path.exists(outdir + dirname):
    OS.makedirs(outdir + dirname)

np.savez(op_fname, opt, snr)
