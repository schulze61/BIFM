#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 12:40:44 2022

@author: joeschulze
"""


import os
import sys

# hack to allow scripts to be placed in subdirectories next to exoplex:
import numpy as np

if not os.path.exists('ExoPlex') and os.path.exists('../ExoPlex'):
    sys.path.insert(1, os.path.abspath('..'))
    
    
import ExoPlex as exo
import ExoPlex.functions as functions

#function to call exoplex for bayesian models

def call_ExoPlex(Mc, FeMg, SiMg, FeO = 0):
    Pressure_range_mantle_UM = '1000 1400000'
    Temperature_range_mantle_UM = '1400 3000'

    Pressure_range_mantle_LM = '1000000 10000000'
    Temperature_range_mantle_LM = '2200 9000'

    core_rad_frac_guess = 0.5
    water_rad_frac_guess = 0.

    combine_phases = True
    use_grids = True
    verbose = False
    
    #Fix minor mantle elements
    CaMg = 0.07
    AlMg = 0.09


    #Other fixed paramets
    wt_frac_water = 0.0
    number_h2o_layers = 0
    water_potential_temp = 300.

    #What fraction of the mantle would you like to be made of FeO? This Fe will be pulled from the core.
    wt_frac_FeO_wanted = FeO #by mass
    Conserve_oxy = False

    #Now we can mix various elements into the core or mantle
    wt_frac_Si_core = 0.0 #by mass <1, note if you conserve oxygen this is calculated for you
    wt_frac_O_core = 0.0 #by mass
    wt_frac_S_core = 0.0 #by mass

    #What potential temperature (in K) do you want to start your mantle adiabat?
    Mantle_potential_temp = 1600.

    #Input the resolution of your upper mantle and lower mantle composition, density grids
    #These are input as number of T, P points. 50 50 = 2500 grid points, which takes about
    #5 minutes to calculate. Lower mantle resolution does not need to be higher since it's
    #mostly ppv.
    resolution_UM = '50 50'
    resolution_LM = '20 20'

    #lastly we need to decide how many layers to put in the planet. This is the resolution of
    #the mass-radius sampling.
    num_mantle_layers = 1000
    num_core_layers = 2000
    
    
    compositional_params = [wt_frac_water,FeMg,SiMg,CaMg,AlMg,wt_frac_FeO_wanted,wt_frac_Si_core, \
                          wt_frac_O_core,wt_frac_S_core,combine_phases,use_grids,Conserve_oxy]
        
        
    #need to calculate mass of planet from Mc and compositional parameters
    cmf_var = functions.get_percents(compositional_params,verbose)[-1]
    Mass_planet = Mc/cmf_var
    
    #print('Masssss: ', Mass_planet)
    #print('Calc CMF: ', cmf_var)
    #print('Fe/Mg  ',  FeMg)
    #print('SiMg  ', SiMg )
    
    
    

    if use_grids == True:
        filename = exo.functions.find_filename(compositional_params,verbose)
    else:
        filename=''

    structure_params =  [Pressure_range_mantle_UM,Temperature_range_mantle_UM,resolution_UM,
                         Pressure_range_mantle_LM, Temperature_range_mantle_LM, resolution_LM,
                         core_rad_frac_guess,Mantle_potential_temp,water_rad_frac_guess,water_potential_temp]


    layers = [num_mantle_layers,num_core_layers,number_h2o_layers]

    #This is where we actually run the planet. First PerPlex grids of mineralogy, density,
    #Cp and alpha are calculated and stored in the Solutions folder. If the file already exists
    #(in name, not necessarily in composition), then PerPlex is not run again.

    Planet = exo.run_planet_mass(Mass_planet,compositional_params,structure_params,layers,filename, verbose)
    
    Mout = (Planet['mass'][-1]/5.97e24)
    Rout = (Planet['radius'][-1]/6371e3)
    
    return Mout, Rout
    
    
    
    
    
    
    
    
    
    
    