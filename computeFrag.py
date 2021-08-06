# -*- coding: utf-8 -*-
"""
K. Trevlopoulos
Based on the MATLAB code by H.Ebrahimian
July 2021

This script reads empirical data and uses them to compute fragility curves
based on MATLAB functions
"""

import matlab.engine
import matlab
import numpy as np
import pandas as pd

# Start the MATLAB engine
eng = matlab.engine.start_matlab()



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# General input

# Enter the filename containing the data
inputFileName = ['IM_and_DS.csv']

df = pd.read_csv(inputFileName[0], sep = ",", engine = 'python')
damageData = np.zeros_like( df )
damageData[:,:] = df
damage_data = matlab.double(damageData.tolist())

# Enter the number of damage scales
NDS = 5
D = matlab.double( range(0,NDS+1) )

# D = {'None','Light','Minor','Moderate','Severe','Collapse'};
# D0 None     -- None
# D1 Light    -- Non-structural damage only
# D2 Minor    -- Significant non-structural damage, minor structural damage
# D3 Moderate -- Significant structural and non-structural damage
# D4 Severe   -- Irreparable structural damage, will require demolition
# D5 Collapse -- Complete structural collapse

# Enter the step and the maximum for the IM vector
dIM = 0.01
IM_max = 10
# The IM vector
IMs = np.arange(dIM,IM_max,dIM)
IM = matlab.double( IMs.tolist() )

# Enther the increment of the vector of the two logistic regression parameters 
dvec_alpha0 = 0.01
dvec_alpha1 = 0.01

# Enter the confidence for the RF
confid = 2
confidence = matlab.double([confid])



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Creating the fragilities Based on logit model

rfragility, sfragility,\
    etaIMc, betaIMc,\
    sample_theta_model =\
    eng.RF_y_logReg_model(damage_data, NDS, D, IM,\
    dvec_alpha0, dvec_alpha1, confidence, nargout=5)

# Conversion back to numpy arrays
rfragility = np.array(rfragility)
sfragility = np.array(sfragility)
etaIMc = np.array(etaIMc)
betaIMc = np.array(betaIMc)
sample_theta_model = np.array(sample_theta_model)
