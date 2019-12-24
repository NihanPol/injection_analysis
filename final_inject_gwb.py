from __future__ import division, print_function

import numpy as np
import sys,os,glob
from collections import OrderedDict
import libstempo as T2
import libstempo.toasim as LT
import argparse
import pickle

#Initialize the parser to accept user inputs
parser = argparse.ArgumentParser(description = "Initiate data sets with stochastic GWB injections with given numpy seed (realization) and amplitudes")

#Required arguments:
parser.add_argument("-seed", nargs = '+', type = int, dest = 'seed', required = True, help = "Seed corresponding to the realization. Can accept multiple values here!")
parser.add_argument("-amps_path", dest = 'amps_path', required = True, help = "Path to numpy file containing array of injected stochastic GWB amplitudes to inject in this dataset; Default: ./injected_amps.npy", default = './injected_amps.npy')
parser.add_argument("-outdir", required = True, help = "Base directory to store output chain and parameter files")

#Optional arguments:
parser.add_argument("--read_pickle", dest = 'read_pickle', action = 'store_true', default = True, help = "Flag to read dataset as pickle to save i/o time; Default: True")
parser.add_argument("--timpath", help = "Path to base directory holding timfiles; Default: /users/nspol/stochastic_11yr_analysis/data/partim/", default = "/users/nspol/stochastic_11yr_analysis/data/partim/")
parser.add_argument("--parpath", help = "Path to directory containing all parfiles; Default: /users/nspol/stochastic_11yr_analysis/data/partim/", default = "/users/nspol/stochastic_11yr_analysis/data/partim/")
parser.add_argument("--ephemeris", dest = 'ephem', help = "Choose solar system ephemeris for loading in pulsars; Default: DE436", choices = ['DE430', 'DE435', 'DE436', 'BayesEphem'], default = 'DE436')
parser.add_argument("--gamma", dest = 'gamma',  help = "Specify index of stochastic GWB powerlaw function; Default: 13./3.", type = float, default = 13./3.)
parser.add_argument("--save_pickle", dest = 'save_pickle', action = 'store_true', default = False, help = "Flag to save dataset as pickle to save i/o time; Default: True")

#Load the arguments:
args = parser.parse_args()

A_gwb = np.load(args.amps_path)

#Need to make this more elegant
if args.read_pickle:
    with open(glob.glob(args.timpath + "*.pkl")[0], 'rb') as f:
        t2psr_orig = pickle.load(f) 

else:
    parfiles = sorted(glob.glob(args.parpath+'*.par'))
    timfiles = sorted(glob.glob(args.timpath+'*.tim'))

    print("Loading in {} and {} par and tim files".format(len(parfiles),len(timfiles)))

#Specify parent directory to hold all the realizations of the different injected amplitudes:
outdir = args.outdir

#This is the meat of the code.
seeds = args.seed

for ii in seeds:
    
    if ii < 0:
        print("Cannot work with seed < 0. Skipping this seed and moving on!")
        continue

    dirname = "/realization_" + str(ii) + '/'
    
    if not os.path.exists(args.outdir + dirname):
        os.makedirs(args.outdir + dirname)
    
    for loc in range(len(A_gwb)):

        #Initialize list to hold pulsar objects
        if args.read_pickle:
            t2psr = list(t2psr_orig)
        else:
            t2psr = []

        if not args.read_pickle:
            for ii in range(len(parfiles)):
                #Read in pulsar parfile and tim file and append to t2psr list
                t2psr.append( T2.tempopulsar(parfile = parfiles[ii], timfile = timfiles[ii],
                                     maxobs=30000, ephem='DE436') )

                #A fail-safe to catch any NaNs or infinites.
                if np.any(np.isfinite(t2psr[ii].residuals())==False)==True:
                    t2psr[ii] = T2.tempopulsar(parfile = parfiles[ii], timfile = timfiles[ii])

                print('\r{0} of {1}'.format(ii+1,len(parfiles)))

        #Inject the GWB into the pulsars:
        LT.createGWB(t2psr, Amp = A_gwb[loc], gam = args.gamma, seed = int(ii))

        amp_dir_name = '/injecting_' + str(loc) + '_gwb/'

        if not os.path.exists(args.outdir + dirname + amp_dir_name):
            os.makedirs(args.outdir + dirname + amp_dir_name)

        if args.save_pickle:
            with open(args.outdir + dirname + amp_dir_name + "inj_psr_pickle.pkl", 'wb') as f:
                pickle.dump(t2psr, f)
        else:
            #save the injected tim files only (not the par files)
            for ii,p in enumerate(t2psr):
                timfile = args.outdir + dirname + amp_dir_name +'{0}.tim'.format(p.name)
                p.savetim(timfile)
