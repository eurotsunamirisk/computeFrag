# computeFrag

https://doi.org/10.5281/zenodo.5167276

A code for computing empirical fragility curves based on logistic regression



This code provides an ensemble of the fragility curves and their corresponding confidence bands for a set of mutually exclusive and collectively exhaustive (MECE) damage states (DS), where DSi: i=0,NDS. The fragility is defined as the probability of damage D exceeding the threshold Di for damage state DSi and is denoted as P(D>Di|IM).  The set of damage levels (Di, i=0:NDS) mark the thresholds of damage states (DSi). 

The parameters of the empirical fragilities associated with different damage level are estimated jointly using Bayesian inference by employing a Markov Chain Monte Carlo Simulation (MCMC) scheme. 



The inputs of the code are the following:

inputFileName: The filename of a csv file containing two columns, one for the intensity measure and one for the damage state. In this example, it is the file IM_and_DS.csv.

NDS: The number of the damage states.

dIM, IM_max: The step and the maximum absolute value for the IM vector.

dvec_alpha0, dvec_alpha1: The increments of the vectors of the two logistic regression parameters.

confid: The number of standard deviations defining the confidence interval for the robust fragility (fixed for all the DS).



The outputs of the code are the following:

sample_theta_model: Ns=500 samples generated from the posterior joint probability distribution for logistic fragility model parameters (e.g., a total of 10 parameters for 5 DS, 2 parameters/DS) using the MCMC procedure. 

rfragility: Robust (the mean fragility among the Ns=500 realizations for each DS) fragility curve with the row showing the IM vector and a number of columns equal to the number of DS.

sfragility: The standard deviation of the Ns=500 fragility curve samples (the same structure as rfragility). For having the fragility with a confidence, you should do this operation: rfragility+confid*sfragility (e.g., if confid=±1, we have 16th and 84th percentiles, if confid=±2, we have 2nd and 98th percentiles).

etaIMc: The median of rfragility (i.e., the IM value corresponding to 50% probability from rfragility). It is a vector with its length equal to the number of DS.

betaIMc: The equivalent logarithmic standard deviation of the rfragility (i.e., half of the logarithmic distance between IM values corresponding to 84% and 16% probabilities, respectively). It is a vector with its length equal to the number of DS.



To execute the code, run the script computeFrag.py.

Requirements: MATLAB, numpy, pandas
