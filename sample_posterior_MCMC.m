%%% MCMC procedure 
%%% Written By: Hossein Ebrahimian, University of Naples Federico II (UNINA)
%%% Laters Updated: 09/2022

function [THETA,accept] = sample_posterior_MCMC(THETA,rank,proposalPDF,proposalPDF_par,priorPDF,priorPDF_par,likelihoodFunction,data)

theta = THETA(rank);

%% sampling new_theta from proposa PDF q(theta) / calculate Proposal Ratio q(old_theta)/q(new_theta)

if strcmp(proposalPDF,'normal')
    
    sigma_theta    = abs(theta)*proposalPDF_par(1);
    new_theta      = normrnd(theta,sigma_theta);  
    proposal_ratio = 1.0;
    
elseif strcmp(proposalPDF,'lognormal')
    
    beta_theta     = proposalPDF_par(1);
    new_theta      = lognrnd(log(theta),beta_theta);   
    proposal_ratio = lognpdf(theta, log(new_theta), beta_theta)/lognpdf(new_theta, log(theta), beta_theta);
    
elseif strcmp(proposalPDF,'uniform')
    
    thetamin  = proposalPDF_par(1);
    thetamax  = proposalPDF_par(2);
    new_theta = unifrnd(thetamin,thetamax); 
    proposal_ratio = 1.0;
    
elseif strcmp(proposalPDF,'kernel')     
    
    seeds     = proposalPDF_par{1,1}(rank,:);
    weights   = proposalPDF_par{1,2};
    new_theta = sampleTheta_kernel(seeds,weights);
    proposal_ratio = calculateKernel(theta,seeds,weights)/calculateKernel(new_theta,seeds,weights);
    
end

%% calculate the ratio of priors

if strcmp(priorPDF,'normal')
    
    meanpriorPDF   = priorPDF_par(1);
    sigmapriorPDF  = priorPDF_par(2);
    ratioPrior     = normpdf(new_theta, meanpriorPDF, sigmapriorPDF)/normpdf(theta, meanpriorPDF, sigmapriorPDF);
    
elseif strcmp(priorPDF,'lognormal')
    
    medianpriorPDF = priorPDF_par(1);
    betapriorPDF   = priorPDF_par(2);
    ratioPrior     = lognpdf(new_theta, log(medianpriorPDF), betapriorPDF)/lognpdf(theta, log(medianpriorPDF), betapriorPDF); 
    
elseif strcmp(priorPDF,'uniform')
    
    ratioPrior     = 1.0;

end

%% calculate the likelihood function 

THETA(rank) = new_theta;
[Likelihoodnew,logLikelihoodnew] = likelihoodFunction(data(1:end-1),THETA,data(end));
THETA(rank) = theta;
[Likelihood,logLikelihood]       = likelihoodFunction(data(1:end-1),THETA,data(end));

%% calculate the acceptance probability

if (isnan(Likelihoodnew/Likelihood) || Likelihoodnew/Likelihood==0)
    pratio =  exp(logLikelihoodnew-logLikelihood)*ratioPrior;
else
    pratio =  (Likelihoodnew/Likelihood)*ratioPrior;
end

alpha = min([1  pratio*proposal_ratio]); 
u = rand;                                    
if u < alpha
    theta = new_theta; 
    accept = 1;
else
    accept = 0;
end
        
THETA(rank) = theta;

%% END
        