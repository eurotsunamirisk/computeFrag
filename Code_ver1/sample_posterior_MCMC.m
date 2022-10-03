%% MCMC procedure 

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
    
    psample   = proposalPDF_par;
    new_theta = psample.random;
    proposal_ratio = ksdensity(psample.InputData.data,theta,'function','pdf')/ksdensity(psample.InputData.data,new_theta,'function','pdf');

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
    
elseif strcmp(priorPDF,'kernel')
  
    ratioPrior     = ksdensity(priorPDF_par.InputData.data,new_theta,'function','pdf')/ksdensity(priorPDF_par.InputData.data,theta,'function','pdf');
    
end

%% calculate the likelihood function 

THETA(rank) = new_theta;
%Likelihoodnew = likelihoodFunction(data,THETA);
[Likelihoodnew,logLikelihoodnew] = likelihoodFunction(data,THETA);
THETA(rank) = theta;
%Likelihood    = likelihoodFunction(data,THETA);
[Likelihood,logLikelihood]      = likelihoodFunction(data,THETA);

%% calculate the acceptance probability

%pratio =  (Likelihoodnew/Likelihood)*ratioPrior;

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
        