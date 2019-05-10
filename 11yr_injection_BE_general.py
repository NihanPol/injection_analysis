from __future__ import division

import numpy as np
import glob, os, json, pickle, sys
import matplotlib.pyplot as plt
import scipy.linalg as sl
import argparse

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

from enterprise_extensions import models, model_utils

#Initialize the parser to accept user inputs
parser = argparse.ArgumentParser(description = "Initiate detection run with model 2A")

#Required arguments:
parser.add_argument("-realiz", required = True, help = "Seed corresponding to the realization")
parser.add_argument("-amp_index", required = True, type = int, help = "Index of the amplitude injected into the data set you're working with")
parser.add_argument("-outdir", required = True, help = "Base directory to store output chain and parameter files")

#Optional arguments:
parser.add_argument("--timpath", help = "Path to base directory holding timfiles of all realizations and amplitude injections; Default: /gpfs/scratch/nspol/real_injected_timfiles/", default = "/gpfs/scratch/nspol/real_injected_timfiles/")
parser.add_argument("--parpath", help = "Path to directory containing all parfiles; Default: /gpfs/home/nspol/stochastic_11yr_analysis/data/partim/", default = "/gpfs/home/nspol/stochastic_11yr_analysis/data/partim/")
parser.add_argument("--noisepath", help = "Path to directory holding noisefiles; Default: /gpfs/home/nspol/stochastic_11yr_analysis/data/noisefiles/", default = "/gpfs/home/nspol/stochastic_11yr_analysis/data/noisefiles/")
parser.add_argument("--ephemeris", dest = 'ephem', help = "Choose solar system ephemeris for loading in pulsars; Default: DE436", choices = ['DE430', 'DE435', 'DE436', 'BayesEphem'], default = 'DE436')
parser.add_argument("--useBE", dest = 'useBE', action = 'store_true', default = False, help = "Flag to use BayesEphem in detection run; Default: False")
parser.add_argument("-ul", "--upper-limit", dest = 'ul', action = 'store_true',  help = "Perform an upper limit run instead of detection run; Default: False", default = False)
parser.add_argument("--gamma", dest = 'gamma',  help = "Specify index of stochastic GWB powerlaw function; Default: 13./3.", type = float, default = 13./3.)
parser.add_argument("--nsamples", dest = 'nsamples', help = 'Number of samples in output chain; Default: 5e6', type = int, default = 500000)
parser.add_argument("--amps_path", dest = 'amps_path', help = "Path to numpy file containing array of injected stochastic GWB amplitudes; Default: ./injected_amps.npy", default = './injected_amps.npy')
parser.add_argument("--psrlist", dest = 'psrlist', default = '', help = "Provide a text file of pulsar names to use in the detection analysis")

parser.add_argument("--dmvar", dest = 'dmvar', default = False, help = "Flag to enable DM variations; Default: False")
parser.add_argument("--dm_gp", dest = 'dm_gp', default = False, help = "Flag to toggle on DM GP; Default: False")
parser.add_argument("--dm_annual", dest = 'dm_annual', default = False, help = "Flag to toggle on DM GP; Default: False")
parser.add_argument("--dm_chrom", dest = 'dm_chrom', default = False, help = "Flag to toggle on DM GP; Default: False")

#Load the arguments:
args = parser.parse_args()

if args.dm_gp:
     dm_gp = 'gp'

#Get the name of directory corresponding to amplitude of gwb:

#A_gwb_1 = np.linspace(1e-16, 1e-15, 20)
#A_gwb_2 = np.linspace(1.1e-15, 5e-15, 10)

#A_gwb = np.append(A_gwb_1, A_gwb_2)
A_gwb = np.load(args.amps_path)
print args

#injection_dir =	'injections/' +	'injecting_' + str(A_gwb[loc]) + '_gwb/'

#print injection_dir

#realiz = np.int(sys.argv[1])

#datadir = '/gpfs/home/nspol/stochastic_11yr_analysis/data/'
#realiz_dir = '/gpfs/scratch/nspol/real_injected_timfiles/realization_' + str(realiz) + '/'
#injection_dir = 'injecting_' + str(A_gwb[loc]) + '_gwb/'

def get_noise_from_enterprise(noisefile):
    
    with open(noisefile) as f:
        params = json.load(f)
    
    return params

parfiles = sorted(glob.glob(args.parpath + '*.par'))
#timfiles = sorted(glob.glob(realiz_dir + injection_dir + '*.tim'))
timfiles = sorted(glob.glob(args.timpath + '/realization_' + args.realiz + '/injecting_' + str(A_gwb[args.amp_index]) + '_gwb/*.tim'))
noisefiles = sorted(glob.glob(args.noisepath + '*.json'))

if args.psrlist:
    psrlist = np.loadtxt(args.psrlist, dtype = 'string')

psrs = []
for p, t in zip(parfiles, timfiles):
    if args.ephem in ['DE430', 'DE435', 'DE436']:
        psr = Pulsar(p, t, ephem=args.ephem) #Cannot read in pulsars with BayesEphem
    else:
        psr = Pulsar(p, t, ephem = 'DE436')
    #Check if this pulsar has baseline > 3 years
    if psr.name in psrlist:
        psrs.append(psr)
    else:
        continue
    
params = {}
for nfile in noisefiles:
    params.update(get_noise_from_enterprise(nfile))
       

pta = models.model_general(psrs, common_psd = 'powerlaw', noisedict = params, gamma_common = args.gamma,
                      upper_limit = args.ul, bayesephem = args.useBE, dm_var = args.dmvar, dm_type = dm_gp, dm_annual = args.dm_annual, dm_chrom = args.dm_chrom)

outdir = args.outdir + '/realization_' + args.realiz + '/injection_' + str(args.amp_index) + '/'
sampler = model_utils.setup_sampler(pta, resume=True, outdir=outdir)

N = int(args.nsamples) # one mega-sample!
x0 = np.hstack(p.sample() for p in pta.params)
sampler.sample(x0, N, AMweight=25, SCAMweight=40, DEweight=55)
