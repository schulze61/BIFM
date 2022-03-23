#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 16:46:10 2022

@author: joeschulze
"""
import numpy as np
import emcee
import matplotlib.pyplot as plt
import scipy.stats as sp
import pandas as pd

import call_exoplex as ce


#define observable parameters
FeMgs = 0.81; sigFeMgs = 0.25;
SiMgs = 0.95; sigSiMgs = 0.25;
Mp = 1.0; sigMp = 0.05
Rp = 1.0; sigRp = 0.01

Mc_guess = 0.33


#define liklihood function
def log_likelihood(theta):
    Mc, FeMg, SiMg = theta
    print(Mc, FeMg, SiMg)
    try:
        model = ce.call_ExoPlex(Mc, FeMg, SiMg)
    except (AssertionError, SystemExit):
        return -np.inf
    
    dM = (model[0]-Mp)
    dR = (model[0]-Rp)
    
    return -0.5 + (dM/sigMp)**2 + (dR/sigRp)**2


#define priors
def prior_FeMg(FeMg):
    a, b = (0.2 - FeMgs) / sigFeMgs, (1000.0 - FeMgs) / sigFeMgs
    #return sp.norm.pdf(FeMg, loc = FeMgs, scale = sigFeMgs)
    return sp.truncnorm.pdf(FeMg, a, b, loc=FeMgs, scale= sigFeMgs)

def prior_SiMg(SiMg):
    a, b = (0.2 - SiMgs) / sigSiMgs, (1000.0 - SiMgs) / sigSiMgs
    #return sp.norm.pdf(FeMg, loc = FeMgs, scale = sigFeMgs)
    return sp.truncnorm.pdf(SiMg, a, b, loc=SiMgs, scale= sigSiMgs)

def prior_Mc(Mc):
    a, b = (0.1 - Mc_guess) / sigMp, (1000.0 - Mc_guess) / sigMp
    #return sp.norm.pdf(FeMg, loc = FeMgs, scale = sigFeMgs)
    #return sp.uniform.pdf(Mc, loc = 0, scale = Mp + 3.0*sigMp)
    return sp.truncnorm.pdf(Mc, a, b, loc=Mc_guess, scale= sigMp)



def log_prior(theta):
    Mc, FeMg, SiMg = theta
    return np.log(prior_FeMg(FeMg)) + np.log(prior_SiMg(SiMg)) + np.log(prior_Mc(Mc))



#now need to define log_prob which is log(priors*likelihood)
def log_prob(theta):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta)


    


nwalkers = 6
Mc_init = sp.uniform.rvs(loc = 0.33*Mp, scale = 0.33*sigMp, size = nwalkers).transpose()
FeMg_init = sp.norm.rvs(loc = FeMgs, scale = 0.1*sigFeMgs, size = nwalkers).transpose()
SiMg_init = sp.norm.rvs(loc = SiMgs, scale = 0.1*sigSiMgs, size = nwalkers).transpose()
pos = np.array([Mc_init, FeMg_init, SiMg_init]).T

ndim = pos.shape[1]


nsteps = 10



sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob, threads = 15)
sampler.run_mcmc(pos, nsteps, progress=True, store = True)

labels = ["Mc", "Fe/Mg", "Si/Mg"]

flat_samples = sampler.get_chain(flat=True)
df = pd.DataFrame(flat_samples, columns  = labels)
df.to_csv('test_mcmc_run_Earth.csv')

