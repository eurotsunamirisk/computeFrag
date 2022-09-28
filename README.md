Guide to accessing the empirical fragility calculation Code

A code for computing empirical fragility curves based on generalized linear regression models. This code provides an ensemble of the fragility curves and their corresponding confidence bands for a set of mutually exclusive and collectively exhaustive (MECE) damage states. To use this code, you should cite the following manuscript:

Jalayer, F., Ebrahimian, H., Trevlopoulos, K., and Bradley, B.: Empirical tsunami fragility modelling for hierarchical damage levels. EGUsphere, https://doi.org/10.5194/egusphere-2022-206, 2022.

The parameters of the empirical fragilities associated with different damage level are estimated jointly using Bayesian inference by employing a Markov Chain Monte Carlo Simulation (MCMC) scheme.

The inputs of the code that should be entered in the main file “main_script_Bayesian_fragility_model.m” are the following:

Filename: The filename of a csv file containing two columns, one for the intensity measure and one for the damage state. 

output_folder: The folder that the outputs will be stored.

excel_filename: The name of the Excel file name for storing the fragility data. 

D: Definition of the damage levels (e.g. 0:2 showing that we have 3 damage levels 0, 1 and 2).

dvec_alpha0, dvec_alpha1: The increments of the vectors of the two logistic regression parameters.

do_MCMC_M1:   Perform MCMC for fragility model 1 (for the models, please see the corresponding paper) (if =0, do not perform MCMC and read MCMC data from output_folder; if =1, do MCMC)

do_MCMC_M2:   Perform MCMC for fragility model 2 (for the models, please see the corresponding paper) (if =0, do not perform MCMC and read MCMC data from output_folder; if =1, do MCMC)

do_MCMC_M3:   Perform MCMC for fragility model 3 (for the models, please see the corresponding paper) (if =0, do not perform MCMC and read MCMC data from output_folder; if =1, do MCMC)

COVprior_M1: Coefficient of variation for the prior parameters of model 1 (please see the corresponding paper) 

COVprior_M2: Coefficient of variation for the prior parameters of model 2 (please see the corresponding paper) 

COVprior_M3: Coefficient of variation for the prior parameters of model 3 (please see the corresponding paper) 

dIM, IM_max: The step and the maximum absolute value for the IM vector.

confidence: The number of standard deviations defining the confidence interval for the robust fragility (fixed for all the DS).

do_modelClassSelection: if =0, do not perform model class selection; if =1, do model class selection (please see the corresponding paper) 

The outputs of the code are stored in the excel files defined by the name excel_filename. 

The output format is described in the file ReadMe file for fragility data.pdf in
https://github.com/eurotsunamirisk/etris_data_and_data_products/tree/main/etris_data_products/Fragility_Curves 
