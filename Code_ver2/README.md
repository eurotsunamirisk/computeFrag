# ComputeFrag V2.0

ComputeFrag V2.0 is a code for computing robust empirical fragility curves using a generalized regression model considering Hierarchical fragility modeling. This code calculates fragility curves based on maximum likelihood estimation methods (basic and hierarchical modeling) or using Bayesian model class selection (BMCS) to estimate fragility curves with their corresponding confidence bands for a set of mutually exclusive and collectively exhaustive damage states and different classes of buildings or infrastructure. The code utilizes Bayesian model class selection (BMCS) to identify the best link model to employ in the generalized linear regression scheme.

**Citation**:
<br> Jalayer, F., Ebrahimian, H., Trevlopoulos, K., and Bradley, B. (2022) Empirical tsunami fragility modelling for hierarchical damage levels. **doi.org/10.5194/egusphere-2022-206.** [Cite](https://doi.org/10.5194/egusphere-2022-206)
<br>
<br> _Please consider citing the following DOI if you directly use ComputeFrag in your work: **doi.org/10.5281/zenodo.5167276.**_ [Cite](https://doi.org/10.5281/zenodo.5167276)

## Files
To run the program use **_main_script_Bayesian_fragility_model.m_** 

## Inputs
The inputs of the code are the following:

- **Input Filename:** The filename of a CSV file containing two columns, one for the intensity measure and one for the damage state.
- **output_folder:** The folder where the outputs will be stored.
- **excel_filename:** The name of the Excel file for storing the fragility data.
- **D:** Definition of the damage levels (e.g., **0:2 showing 3 damage levels 0, 1, and 2**).
- **dvec_alpha0, dvec_alpha1:** The increments of the vectors of the two logistic regression parameters.
- **do_MCMC_M1, do_MCMC_M2, do_MCMC_M3:** Perform MCMC for fragility models 1, 2, and 3 (if 0, do not perform MCMC; if 1, perform MCMC).
- **COVprior_M1, COVprior_M2, COVprior_M3:** Coefficient of variation for the prior parameters of models 1, 2, and 3.
- **dIM, IM_max:** The step and the maximum absolute value for the IM vector.
- **confidence:** The number of standard deviations defining the confidence interval for the robust fragility (fixed for all the DS).
- **do_modelClassSelection:** If 0, do not perform model class selection; if 1, do model class selection.

## Outputs
The output data is available at _.csv_ or _.xlsx_ format, including: 

- **flow depth [m]:** The tsunami fragility curve employs the tsunami flow depth in meters as the measure of intensity.
- **mean-1sigma fragility:** The expected value of the fragility over the vector of fragility parameters for a given intensity minus 1 standard deviation (the lower confidence interval of the fragility).
- **mean fragility:** The expected value of the fragility over the vector of fragility parameters for a given intensity.
- **mean+1sigma fragility:** The expected value of the fragility over the vector of fragility parameters for a given intensity plus 1 standard deviation (the upper confidence interval of the fragility).
- **median:** The median (50%) of the mean fragility curve.
- **logarithmic standard deviation:** The logarithmic standard deviation (dispersion) of the mean fragility curve.
- **epistemic uncertainty:** Epistemic uncertainty due to the uncertainty in the fragility model parameters on the median fragility.

**All of these columns are repeated for different damage scales.**

## Docker repository

<p style="text-align: justify;font-size: 17px" class="has-poppins-font-family">
i) Pull the image from: <a href="https://hub.docker.com/r/eurotsunamirisk/bayesian-fragility-standalone-app"> Docker hub</a> 
<br>ii) Create a folder for the input and the outputs of the code. (e.g. C:\input_output\ (Windows) or /home/user/input_output (Linux)).
<br>iii) Place the input file (e.g. the file building_class_1.csv from the repository). 
<br>
iv) Execute the application using one of the following commands, depending on your operating system:
</p>
<p style="font-size: 15px" class="has-poppins-font-family">
<code>docker run --rm -v C:\input_output:/tmp bayesian-fragility-standalone-app/tmp/building_class_1.csv</code><br><code>docker run --rm -v /home/user/input_output:/tmp bayesian-fragility-standalone-app/tmp/building_class_1.csv</code>
</p>
