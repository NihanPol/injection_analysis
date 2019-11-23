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
import sample_utils as su

from enterprise_extensions import models, model_utils

def get_noise_from_enterprise(noisefile):
    
    with open(noisefile) as f:
        params = json.load(f)
    
    return params

#Initialize the parser to accept user inputs
parser = argparse.ArgumentParser(description = "Single pulsar noise runs with model_singlepsr_noise")

#Required arguments:
parser.add_argument("-realiz", required = True, help = "Seed corresponding to the realization")
parser.add_argument("-amp_index", required = True, type = int, help = "Index of the amplitude injected into the data set you're working with")
parser.add_argument("-psr_index", required = True, type = int, help = "Index of the n'th psr in PTA")
parser.add_argument("-outdir", required = True, help = "Base directory to store output chain and parameter files")

#Optional arguments:
parser.add_argument("--timpath", help = "Path to base directory holding timfiles of all realizations and amplitude injections; Default: /gpfs/scratch/nspol/real_injected_timfiles/", default = "/gpfs/scratch/nspol/real_injected_timfiles/")
parser.add_argument("--parpath", help = "Path to directory containing all parfiles; Default: /gpfs/home/nspol/stochastic_11yr_analysis/data/partim/", default = "/gpfs/home/nspol/stochastic_11yr_analysis/data/partim/")
#parser.add_argument("--noisepath", help = "Path to directory holding noisefiles; Default: /gpfs/home/nspol/stochastic_11yr_analysis/data/noisefiles/", default = "/gpfs/home/nspol/stochastic_11yr_analysis/data/noisefiles/")
parser.add_argument("--ephemeris", dest = 'ephem', help = "Choose solar system ephemeris for loading in pulsars; Default: DE436", choices = ['DE430', 'DE435', 'DE436', 'BayesEphem'], default = 'DE436')
#parser.add_argument("--useBE", dest = 'useBE', action = 'store_true', default = False, help = "Flag to use BayesEphem in detection run; Default: False")
parser.add_argument("-ul", "--upper-limit", dest = 'ul', action = 'store_true',  help = "Perform an upper limit run instead of detection run, with ALL priors set to UL priors; Default: False", default = False)
parser.add_argument("--gamma", dest = 'gamma',  help = "Specify index of stochastic GWB powerlaw function; Default: None", default = 'None')
parser.add_argument("--nsamples", dest = 'nsamples', help = 'Number of samples in output chain; Default: 5e6', type = int, default = 500000)
parser.add_argument("--thin", dest = 'thin', help = 'thinning factor (keep every [thin]th sample)', type = int, default = 10)
parser.add_argument("--amps_path", dest = 'amps_path', help = "Path to numpy file containing array of injected stochastic GWB amplitudes; Default: ./injected_amps.npy", default = './injected_amps.npy')
parser.add_argument("--psrlist", dest = 'psrlist', default = '', help = "Provide a text file of pulsar names to use in the detection analysis")
#parser.add_argument("--orf", dest = 'orf', default = None, help = "Set the orf from None(no spat. corr), monopole, dipole, hd; Default: None")

#Set type of spectrum:
#parser.add_argument("--common_psd", dest = 'common_psd', default = 'powerlaw', help = "Set common red noise spectrum type (powerlaw, spectrum, turnover); Default: powerlaw")
parser.add_argument("--red_psd", dest = 'red_psd', default = 'powerlaw', help = "Set ind. red noise spectrum type (powerlaw, spectrum, turnover); Default: powerlaw")

parser.add_argument("--dm_var", dest = 'dm_var', action = 'store_true', default = False, help = "Flag to enable DM variations; Default: False")
parser.add_argument("--dm_gp", dest = 'dm_gp', action = 'store_true', default = False, help = "Flag to toggle on DM GP; Default: False")
parser.add_argument("--dm_annual", dest = 'dm_annual', action = 'store_true', default = False, help = "Flag to toggle on DM GP; Default: False")
parser.add_argument("--dm_chrom", dest = 'dm_chrom', action = 'store_true', default = False, help = "Flag to toggle on DM GP; Default: False")

#Load the arguments:
args = parser.parse_args()

if args.dm_gp:
    dm_gp = 'gp'

if args.gamma == 'None' or args.gamma is None:
    gamma = None
else:
    gamma = float(args.gamma)

if args.psrlist:
    psrlist = np.loadtxt(args.psrlist, dtype = 'string')

#Get psr_index:
psr_index = int(args.psr_index)

#Pick out pulsar corresponding to psr_index from psrlist:
psr_choice = psrlist[psr_index]

#Get the name of directory corresponding to amplitude of gwb:
A_gwb = np.load(args.amps_path)

parfiles = sorted(glob.glob(args.parpath + '*' + psr_choice + '*.par'))[0]
timfiles = sorted(glob.glob(args.timpath + '/realization_' + args.realiz + '/injecting_' + str(A_gwb[args.amp_index]) + '_gwb/' + '*' + psr_choice + '*.tim'))[0]
#noisefiles = sorted(glob.glob(args.noisepath + '*.json'))

#Read in the par and tim files and create the list of pulsars
psrs = []

if args.ephem in ['DE430', 'DE435', 'DE436']:
    psr = Pulsar(parfiles, timfiles, ephem=args.ephem) #Cannot read in pulsars with BayesEphem
else:
    psr = Pulsar(p, t, ephem = 'DE436')
       
#Setup the PTA
pta = models.model_singlepsr_noise(psrs, tm_var = False, gamma_val = gamma, red_var = True, white_vary = True, 
                                   psd = args.red_psd, upper_limit = args.ul, 
                                   dm_var = args.dm_var, dm_type = dm_gp, dm_annual = args.dm_annual, dm_chrom = args.dm_chrom,)

#Set output directory
outdir = args.outdir + '/realization_' + args.realiz + '/injection_' + str(args.amp_index) + '/' + psr_choice + '/'

sampler = model_utils.setup_sampler(pta, resume=True, outdir=outdir)

N = int(args.nsamples) # one mega-sample!
x0 = np.hstack(p.sample() for p in pta.params)
sampler.sample(x0, N, AMweight=25, SCAMweight=40, DEweight=55)