from __future__ import division,print_function

import numpy as np
import sys,os,glob
from collections import OrderedDict
import libstempo as T2
import libstempo.toasim as LT
import argparse

parpath = '/users/nspol/stochastic_11yr_analysis/data/partim/'
timpath = '/users/nspol/stochastic_11yr_analysis/data/partim/'
noisepath = '/users/nspol/stochastic_11yr_analysis/data/nano_11_noisefiles_DE430_30/'

parfiles = sorted(glob.glob(parpath+'*.par'))
timfiles = sorted(glob.glob(timpath+'*.tim'))
noisefiles = sorted(glob.glob(noisepath+'*.txt'))

print(len(parfiles),len(timfiles))

## Pulsar noise info class
class psr_noise(object):

    def __init__(self, t2psr, noisefile):

        self.name = t2psr.name

        with open(noisefile, 'r') as content_file:
            file_content = content_file.read()

        noiselines = file_content.split('\n')
        self.efacs = []
        self.equads = []
        self.ecorrs = []
        for ll in noiselines:
            if 'efac' in ll:
                self.efacs.append([ll.split()[0].split('efac-')[1],
                              np.double(ll.split()[1])])
            if 'equad' in ll:
                self.equads.append([ll.split()[0].split('equad-')[1],
                               10.0**np.double(ll.split()[1])])
            if 'jitter' in ll:
                self.ecorrs.append([ll.split()[0].split('jitter_q-')[1],
                               10.0**np.double(ll.split()[1])])

        self.efacs = OrderedDict(self.efacs)
        self.equads = OrderedDict(self.equads)
        self.ecorrs = OrderedDict(self.ecorrs)

        for ll in noiselines:
            if 'RN-Amplitude' in ll:
                self.Redamp = 10.0**np.double(ll.split()[1])
            if 'RN-spectral-index' in ll:
                self.Redind = np.double(ll.split()[1])

outdir = '/scratch/nspol/simulated_data/injected_timfiles/'

A_gwb_1 = np.linspace(1e-16, 1e-15, 20)
A_gwb_2 = np.linspace(1.1e-15, 5e-15, 10)
A_gwb = np.append(A_gwb_1, A_gwb_2)

seed_efac = 3
seed_equad = 4
seed_jitter = 5
seed_red = 6
seed_gwb = 7

#for jj in range(10):
for jj in [1000, 2000, 3000, 1993]:    
    for loc in range(len(A_gwb)):
        
        t2psr = []
        
        for ii in range(len(parfiles)):

            t2psr.append( T2.tempopulsar(parfile = parfiles[ii], timfile = timfiles[ii],
                                         maxobs=30000, ephem='DE436') )

            if np.any(np.isfinite(t2psr[ii].residuals())==False)==True:
                t2psr[ii] = T2.tempopulsar(parfile = parfiles[ii], timfile = timfiles[ii])

            print('\r{0} of {1}'.format(ii+1,len(parfiles)))

        pnoise = OrderedDict({})
        for ii, p in enumerate(t2psr):
            pnoise[p.name] = psr_noise(p,noisefiles[ii])

        for ii,p in enumerate(t2psr):

            ## make ideal
            LT.make_ideal(p)

            ## add efacs
            LT.add_efac(p, efac =list(pnoise[p.name].efacs.values()),
                        flagid = 'f', flags = list(pnoise[p.name].efacs.keys()),
                        seed = seed_efac + ii)

            ## add equads
            LT.add_equad(p, equad = list(pnoise[p.name].equads.values()),
                         flagid = 'f', flags = list(pnoise[p.name].equads.keys()),
                         seed = seed_equad + ii)

            ## add jitter
            LT.add_jitter(p, ecorr = list(pnoise[p.name].ecorrs.values()),
                          flagid='f', flags = list(pnoise[p.name].ecorrs.keys()),
                          coarsegrain = 1.0/86400.0, seed=seed_jitter + ii)

            ## add red noise
            LT.add_rednoise(p, pnoise[p.name].Redamp, pnoise[p.name].Redind,
                            components = 30, seed = seed_red + ii)

            #Try fitting the residuals as a test for source of excess signal. Comment out in general.
            #p.fit(iters = 5)

            print(ii+1, p.name, ' ')
            
        print("Injecting amplitude {} with seed {}".format(A_gwb[loc], jj))
        LT.createGWB(t2psr, Amp = A_gwb[loc], gam=13./3., seed = int(jj))

        #outdir = '/users/nspol/stochastic_11yr_analysis/data/injections/'
        dirname = "realization_" + str(jj) + '/'
        amp_dir_name = 'injecting_' + str(A_gwb[loc]) + '_gwb/'

        if not os.path.exists(outdir + dirname):
            os.makedirs(outdir + dirname)

        if not os.path.exists(outdir + dirname + amp_dir_name):
            os.makedirs(outdir + dirname + amp_dir_name)

        for ii,p in enumerate(t2psr):
            timfile = outdir + dirname + amp_dir_name +'{0}.tim'.format(p.name)
            p.savetim(timfile)
