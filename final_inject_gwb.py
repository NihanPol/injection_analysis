from __future__ import division,print_function

import numpy as np
import sys,os,glob
from collections import OrderedDict
import libstempo as T2
import libstempo.toasim as LT
import argparse

#Specify paths to par and tim files. EDIT
parpath = '/users/nspol/stochastic_11yr_analysis/data/partim/'
timpath = '/users/nspol/stochastic_11yr_analysis/data/partim/'

parfiles = sorted(glob.glob(parpath+'*.par'))
timfiles = sorted(glob.glob(timpath+'*.tim'))

print(len(parfiles),len(timfiles))

#Define different stochastic GWB amplitudes to be injected:
A_gwb_1 = np.linspace(1e-16, 1e-15, 20)
A_gwb_2 = np.linspace(1.1e-15, 5e-15, 10)

A_gwb = np.append(A_gwb_1, A_gwb_2)

#Here make 10 different realizations of the GWB by changing the seed:

#Specify parent directory to hold all the realizations of the different injected amplitudes:
outdir = '/scratch/nspol/real_injected_timfiles/'

#This is the meat of the code.
#The list in the for statement contains pseudo-randomly chosen numpy.seed values
for ii in [0, 1, 2, 3, 4, 5, 1000, 1993, 2000, 3000]:
    
    dirname = "realization_" + str(ii) + '/'
    
    if not os.path.exists(outdir + dirname):
        os.makedirs(outdir + dirname)
    
    for loc in range(len(A_gwb)):

        #Initialize list to hold pulsar objects
        t2psr = []
        
        for ii in range(len(parfiles)):
            #Read in pulsar parfile and tim file and append to t2psr list
            t2psr.append( T2.tempopulsar(parfile = parfiles[ii], timfile = timfiles[ii],
                                     maxobs=30000, ephem='DE436') )

            #A fail-safe to catch any NaNs or infinites.
            if np.any(np.isfinite(t2psr[ii].residuals())==False)==True:
                t2psr[ii] = T2.tempopulsar(parfile = parfiles[ii], timfile = timfiles[ii])

            print('\r{0} of {1}'.format(ii+1,len(parfiles)))

        #Inject the GWB into the pulsars:
        LT.createGWB(t2psr, Amp=A_gwb[loc], gam=13./3., seed = int(ii))

        #outdir = '/users/nspol/stochastic_11yr_analysis/data/injections/'
        amp_dir_name = 'injecting_' + str(A_gwb[loc]) + '_gwb/'

        if not os.path.exists(outdir + dirname + amp_dir_name):
            os.makedirs(outdir + dirname + amp_dir_name)

        #save the injected tim files only (not the par files)
        for ii,p in enumerate(t2psr):
            timfile = outdir + dirname + amp_dir_name +'{0}.tim'.format(p.name)
            p.savetim(timfile)
