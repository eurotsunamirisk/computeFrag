%%% This is the main script for Bayesian Empirical fragility modelling
%%% Written By: Hossein Ebrahimian, University of Naples Federico II (UNINA)
%%% Laters Updated: 09/2022

clc
clear

%% General Input 

%%% -------------- Input/output files
filename  = 'Palu Building Class 1.csv';
output_folder = 'MCMC_output_BC1_2018';
excel_filename = 'xxx';

%%% --------------  Definition of D (Damage levels)
D = 0:2;   

%%% -------------- Input parameters for MCMC
%%% Increment of the vector of the two logistic regression parameter for each DS
dvec_alpha0 = 0.01;
dvec_alpha1 = 0.01;

%%% Do MCMC
do_MCMC_M1   = 0;
COVprior_M1  = 0.80*4*ones(4,1);   

do_MCMC_M2   = 0;
COVprior_M2  = 0.80*4*ones(4,1);   

do_MCMC_M3   = 0;
COVprior_M3  = 0.80*4*[12,1,6,1];    

%%% -------------- Input parameters for Fragility Curve
%%% vector of Sa for Fragility plot
dIM = 0.01;
IM_max = 10;

%%% For RF
confidence = 1;

do_modelClassSelection = 1;

%% Create OUTPUT folder

if exist(output_folder,'dir')==0
    mkdir(output_folder)
end

%% Vector IM

IM = [1e-6,dIM:dIM:IM_max];

%% Damage state Definition for program

NDS=D(end)-D(1);   % (=length(D)-1); NDS is the real numer of (DS)-1 / number of distances between damage states (i.e., number of DS that the model parameters are estimated for) 

%%% For Plot
vecDS = {};
vecDSexcel = {};
for i=2:length(D)
    vecDS = [vecDS,{['\itD\rm_{',num2str(D(i)),'}']}];
    vecDSexcel = [vecDSexcel,{['D',num2str(D(i))]}];
end

%% Reading DS

fid = fopen(filename);
damage_data = textscan(fid,'%f %f','delimiter',',');
fclose(fid);

DS = cell(1,NDS+1);

for j=1:NDS+1
    DS{1,j} = damage_data{1,1}(damage_data{1,2}==D(j));    % DS is a NDS+1=length(D) cell and each cell has the IM data inside 
end

%% Model parameter estimation for Model 1

GRM = 'logit';          

display(['------ Model 1, using ',GRM,' link function ------']) 
display(' ') 

theta_prior_model1 = calculate_glm(DS,NDS,IM,GRM);

if do_MCMC_M1

[sample_theta_model1,priorPDF_par_model1] = posteriorFragilityFunction(DS,dvec_alpha0,dvec_alpha1,theta_prior_model1,'Blockwise',GRM,COVprior_M1);

save([output_folder,'\MCMC_M1'],'sample_theta_model1','priorPDF_par_model1');

else
    
load([output_folder,'\MCMC_M1'])    

end

%% Model parameter estimation for Model 2

GRM = 'probit';           

display(' ') 
display(['------ Model 2, using ',GRM,' link function ------']) 
display(' ') 

theta_prior_model2 = calculate_glm(DS,NDS,IM,GRM);

if do_MCMC_M2

[sample_theta_model2,priorPDF_par_model2] = posteriorFragilityFunction(DS,dvec_alpha0,dvec_alpha1,theta_prior_model2,'Blockwise',GRM,COVprior_M2);

save([output_folder,'\MCMC_M2'],'sample_theta_model2','priorPDF_par_model2');

else
    
load([output_folder,'\MCMC_M2'])

end

%% Model parameter estimation for Model 3

GRM = 'comploglog';          

display(' ') 
display(['------ Model 3, using ',GRM,' link function ------']) 
display(' ') 

theta_prior_model3 = calculate_glm(DS,NDS,IM,GRM);

if do_MCMC_M3

[sample_theta_model3,priorPDF_par_model3] = posteriorFragilityFunction(DS,dvec_alpha0,dvec_alpha1,theta_prior_model3,'Blockwise',GRM,COVprior_M3);

save([output_folder,'\MCMC_M3'],'sample_theta_model3','priorPDF_par_model3');

else
    
load([output_folder,'\MCMC_M3'])

end

%% Estimate the Robust Fragility / FA for Model 1 

GRM = 'logit'; 

tol = -1e-6;
RF_model1 = RobustFragility_DS(sample_theta_model1,IM,NDS,confidence,GRM,tol);

FA_model1 = calculate_fragility(IM',theta_prior_model1,NDS,GRM);
    
%% Estimate the Robust Fragility / FA for Model 2 

GRM = 'probit';          

RF_model2 = RobustFragility_DS(sample_theta_model2,IM,NDS,confidence,GRM,tol);

FA_model2 = calculate_fragility(IM',theta_prior_model2,NDS,GRM);

%% Estimate the Robust Fragility / FA for Model 3 

GRM = 'comploglog';           

RF_model3 = RobustFragility_DS(sample_theta_model3,IM,NDS,confidence,GRM,tol);

FA_model3 = calculate_fragility(IM',theta_prior_model3,NDS,GRM);

%% Creating the fragilities based on Basic Method for Model 1

GRM = 'logit';           

[theta_basic_model1,DSj,vecNcol_NDS] = calculate_glm_basic(DS,NDS,IM,GRM);

FA_basic_model1 = calculate_fragility_basic(IM',theta_basic_model1,NDS,GRM);

%% Creating the fragilities based on Basic Method for Model 1

GRM = 'probit';          

theta_basic_model2 = calculate_glm_basic(DS,NDS,IM,GRM);

FA_basic_model2 = calculate_fragility_basic(IM',theta_basic_model2,NDS,GRM);

%% Creating the fragilities based on Basic Method for Model 1

GRM = 'comploglog';           

theta_basic_model3 = calculate_glm_basic(DS,NDS,IM,GRM);

FA_basic_model3 = calculate_fragility_basic(IM',theta_basic_model3,NDS,GRM);

%% Save

RF2ExcelFile(output_folder,excel_filename,IM',RF_model1,vecDSexcel,'M1',1)

RF2ExcelFile(output_folder,excel_filename,IM',RF_model2,vecDSexcel,'M2',1)

RF2ExcelFile(output_folder,excel_filename,IM',RF_model3,vecDSexcel,'M3',1)

%% Model Class Selection

if do_modelClassSelection 

display(' ') 
display('------ Calculating Evidence for Model 1 ------')
display(' ') 
GRM = 'logit';           
[log_evidence(1),meanlogLikelihood(1),meanlogratioP(1)] = calculate_logE(RF_model1.sample_theta,@likelihoodFunction_DS,DS,priorPDF_par_model1,'normal',GRM);

display(' ') 
display('------ Calculating Evidence for Model 2 ------')
display(' ') 
GRM = 'probit';          
[log_evidence(2),meanlogLikelihood(2),meanlogratioP(2)] = calculate_logE(RF_model2.sample_theta,@likelihoodFunction_DS,DS,priorPDF_par_model2,'normal',GRM);

display(' ') 
display('------ Calculating Evidence for Model 3 ------')
display(' ') 
GRM = 'comploglog';           
[log_evidence(3),meanlogLikelihood(3),meanlogratioP(3)] = calculate_logE(RF_model3.sample_theta,@likelihoodFunction_DS,DS,priorPDF_par_model3,'normal',GRM);

w = [1/3,1/3,1/3];   
P_M = exp(log_evidence).*w;
P_M = P_M/sum(P_M);

display(['------ The posterior probability of Model Class 1 is ',num2str(P_M(1),'%4.3f')])
display(['------ The posterior probability of Model Class 2 is ',num2str(P_M(2),'%4.3f')])
display(['------ The posterior probability of Model Class 3 is ',num2str(P_M(3),'%4.3f')])

end

%% END
