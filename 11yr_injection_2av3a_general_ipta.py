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
parser.add_argument("-ul", "--upper-limit", dest = 'ul', action = 'store_true',  help = "Perform an upper limit run instead of detection run, with ALL priors set to UL priors; Default: False", default = False)
parser.add_argument("--gamma", dest = 'gamma',  help = "Specify index of stochastic GWB powerlaw function; Default: 13./3.", default = 13./3.)
parser.add_argument("--nsamples", dest = 'nsamples', help = 'Number of samples in output chain; Default: 5e6', type = int, default = 500000)
parser.add_argument("--thin", dest = 'thin', help = 'thinning factor (keep every [thin]th sample)', type = int, default = 10)
parser.add_argument("--amps_path", dest = 'amps_path', help = "Path to numpy file containing array of injected stochastic GWB amplitudes; Default: ./injected_amps.npy", default = './injected_amps.npy')
parser.add_argument("--psrlist", dest = 'psrlist', default = '', help = "Provide a text file of pulsar names to use in the detection analysis")

#Set type of spectrum:
parser.add_argument("--common_psd", dest = 'common_psd', default = 'powerlaw', help = "Set common red noise spectrum type (powerlaw, spectrum, turnover); Default: powerlaw")
parser.add_argument("--red_psd", dest = 'red_psd', default = 'powerlaw', help = "Set ind. red noise spectrum type (powerlaw, spectrum, turnover); Default: powerlaw")

#Set priors explicitly here:
parser.add_argument("--upper_limit_red", dest = 'upper_limit_red', action = 'store_true', default = None, help = 'Flag to set ONLY individual red noise prior to UL prior; Default: off')
parser.add_argument("--upper_limit_dm", dest = 'upper_limit_dm', action = 'store_true', default = None, help = 'Flag to set ONLY individual DM GP prior to UL prior; Default: off')
parser.add_argument("--upper_limit_common", dest = 'upper_limit_common', action = 'store_true', default = None, help = 'Flag to set ONLY common red noise prior to UL prior; Default: off')

parser.add_argument("--dm_var", dest = 'dm_var', action = 'store_true', default = False, help = "Flag to enable DM variations; Default: False")
parser.add_argument("--dm_gp", dest = 'dm_gp', action = 'store_true', default = False, help = "Flag to toggle on DM GP; Default: False")
parser.add_argument("--dm_annual", dest = 'dm_annual', action = 'store_true', default = False, help = "Flag to toggle on DM GP; Default: False")
parser.add_argument("--dm_chrom", dest = 'dm_chrom', action = 'store_true', default = False, help = "Flag to toggle on DM GP; Default: False")

#Load the arguments:
args = parser.parse_args()

if args.dm_gp:
    dm_gp = 'gp'

if args.gamma == 'None':
    gamma = None
else:
    gamma = float(args.gamma)

#Get the name of directory corresponding to amplitude of gwb:

A_gwb = np.load(args.amps_path)

parfiles = sorted(glob.glob(args.parpath + '*.par'))
timfiles = sorted(glob.glob(args.timpath + '/realization_' + args.realiz + '/injecting_' + str(args.amp_index) + '_gwb/*.tim'))
noisefiles = sorted(glob.glob(args.noisepath + '*.json'))

if args.psrlist:
    psrlist = np.loadtxt(args.psrlist, dtype = 'string')

#Read in the par and tim files and create the list of pulsars
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

#Read in and store the noise properties
noise_params = {}
for nfile in noisefiles:
    noise_params.update(get_noise_from_enterprise(nfile))


#Initialize dictionary with two keys for model 2a and 3a
sub_pta = dict.fromkeys(np.arange(0, 2))

sub_pta[0] = models.model_general(psrs, noisedict = noise_params, gamma_common = gamma, orf = None, common_psd = args.common_psd, red_psd = args.red_psd,
                      upper_limit = args.ul, bayesephem = args.useBE, dm_var = args.dm_var, dm_type = dm_gp, dm_annual = args.dm_annual, dm_chrom = args.dm_chrom, upper_limit_red = args.upper_limit_red, upper_limit_dm = args.upper_limit_dm, upper_limit_common = args.upper_limit_common)

sub_pta[1] = models.model_general(psrs, noisedict = noise_params, gamma_common = gamma, orf = 'hd', common_psd = args.common_psd, red_psd = args.red_psd,
                      upper_limit = args.ul, bayesephem = args.useBE, dm_var = args.dm_var, dm_type = dm_gp, dm_annual = args.dm_annual, dm_chrom = args.dm_chrom, upper_limit_red = args.upper_limit_red, upper_limit_dm = args.upper_limit_dm, upper_limit_common = args.upper_limit_common)

#This is the supermodel initialized here.
#I call it PTA so as to have to change the jump proposal code below.
pta = model_utils.HyperModel(sub_pta)

#Set output directory
outdir = args.outdir + '/realization_' + args.realiz + '/injection_' + str(args.amp_index) + '/'

sampler = pta.setup_sampler(resume=True, outdir=outdir)

N = int(5e6) # one mega-sample!
x0 = pta.initial_sample()
sampler.sample(x0, N, AMweight=25, SCAMweight=40, DEweight=55)

#This does not play well with the hypermodel setup. 
#Commenting out for now, will come back to it later
#Use the above default sampler setup w/o specific jump proposals
"""
#Setup up groupings for IPTA sampling

#IPTA grouping and sampling:

x0 = np.hstack([noise_params[p.name] if p.name in noise_params.keys()
                else p.sample() for p in pta.params])  # initial point
ndim = len(x0)

# set initial cov stdev to (starting order of magnitude)/10
stdev = np.array([10 ** np.floor(np.log10(abs(x))) / 10 for x in x0])
cov = np.diag(stdev**2)

# generate custom sampling groups
groups = [list(range(ndim))]

# pulsar noise groups (RN + DM)
for psr in psrs:
    this_group = [pta.param_names.index(par)
                  for par in pta.param_names if psr.name in par]
    groups.append(this_group)

groups.append([pta.param_names.index('gw_log10_A')])

if args.useBE:
    # all BE params
    this_group = [pta.param_names.index(par)
                  for par in pta.param_names
                  if 'jup_orb' in par or 'mass' in par or 'frame_drift' in par]
    groups.append(this_group)

    # jup_orb elements + GWB
    this_group = [pta.param_names.index(par)
                  for par in pta.param_names if 'jup_orb' in par]
    this_group.append(pta.param_names.index('gw_log10_A'))
    groups.append(this_group)

sampler = ptmcmc(ndim, pta.get_lnlikelihood, pta.get_lnprior, cov, groups=groups,
                 outDir=outdir, resume=True)

#Write out the pars file for easyness later:
np.savetxt(outdir+'/pars.txt',
               list(map(str, pta.param_names)), fmt='%s')

# additional proposals
full_prior = su.build_prior_draw(pta, pta.param_names, name='full_prior')
sampler.addProposalToCycle(full_prior, 10)

RNA_params = [pname for pname in pta.param_names if 'red_noise_log10_A' in pname]
RN_loguni = su.build_loguni_draw(pta, RNA_params, (-20,-11), name='RN_loguni')
sampler.addProposalToCycle(RN_loguni, 5)

GWB_loguni = su.build_loguni_draw(pta, 'gw_log10_A', (-18,-12), name='GWB_loguni')
sampler.addProposalToCycle(GWB_loguni, 5)

if args.useBE:
    BE_params = [pta.param_names[ii] for ii in groups[6]]
    BE_prior = su.build_prior_draw(pta, BE_params, name='BE_prior')
    sampler.addProposalToCycle(BE_prior, 5)

# SAMPLE!!
Nsamp = args.nsamples * args.thin
sampler.sample(x0, Nsamp,
               SCAMweight=30, AMweight=20, DEweight=50,
burn=int(5e4), thin=args.thin)
"""
