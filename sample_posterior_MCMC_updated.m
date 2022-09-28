%%% MCMC procedure 
%%% Written By: Hossein Ebrahimian, University of Naples Federico II (UNINA)
%%% Laters Updated: 09/2022

function [THETA,accept] = sample_posterior_MCMC_updated(THETA,seeds,weights,priorPDF,priorPDF_par,weights_prior,likelihoodFunction,data)

%% sampling new_theta from proposa PDF q(theta) / calculate Proposal Ratio q(old_theta)/q(new_theta)

new_THETA = sampleTheta_kernel(seeds,weights);

if strcmp(priorPDF{1,1},'lognormal')
    while any(new_THETA<0)
        new_THETA = sampleTheta_kernel(seeds,weights);
    end
end

%proposal_ratio = kernelPDF(THETA,seeds,'adaptive')/kernelPDF(new_THETA,seeds,'adaptive');

proposal_ratio = calculateKernel(THETA,seeds,weights)/calculateKernel(new_THETA,seeds,weights);

%% calculate the ratio of priors

if ~strcmp(priorPDF{1,1},'kernel')

    if ~all(strcmp(priorPDF,'normal') | strcmp(priorPDF,'lognormal'))

        ratioPrior = zeros(1,length(THETA));

        for i=1:length(THETA)

            if strcmp(priorPDF{1,i},'normal')

            meanpriorPDF   = priorPDF_par{1,i}(1);
            sigmapriorPDF  = priorPDF_par{1,i}(2);
            ratioPrior(i)  = normpdf(new_THETA(i), meanpriorPDF, sigmapriorPDF)/normpdf(THETA(i), meanpriorPDF, sigmapriorPDF);

            elseif strcmp(priorPDF{1,i},'lognormal')

            medianpriorPDF = priorPDF_par{1,i}(1);
            betapriorPDF   = priorPDF_par{1,i}(2);
            ratioPrior(i)  = lognpdf(new_THETA(i), log(medianpriorPDF), betapriorPDF)/lognpdf(THETA(i), log(medianpriorPDF), betapriorPDF); 

            elseif strcmp(priorPDF{1,i},'uniform')

            ratioPrior(i)  = 1.0;

            end

        end
        
        ratioPrior = prod(ratioPrior);

    else

        priorNormal_par = cell2mat(priorPDF_par');
        ratioPrior = calculateMVN(new_THETA,priorNormal_par(:,1),priorNormal_par(:,2),eye(length(new_THETA)),priorPDF{1,1})/calculateMVN(THETA,priorNormal_par(:,1),priorNormal_par(:,2),eye(length(new_THETA)),priorPDF{1,1});

    end

else

    samplePrior = cell2mat(priorPDF_par);    
    ratioPrior  = calculateKernel(new_THETA,samplePrior,weights_prior)/calculateKernel(THETA,samplePrior,weights_prior);

end

%% calculate the likelihood function 

[Likelihoodnew,logLikelihoodnew] = likelihoodFunction(data(1:end-1),new_THETA,data(end));
[Likelihood,logLikelihood]       = likelihoodFunction(data(1:end-1),THETA,data(end));

%% calculate the acceptance probability

if (isnan(Likelihoodnew/Likelihood) || Likelihoodnew/Likelihood==0)
    pratio =  exp(logLikelihoodnew-logLikelihood)*ratioPrior;
else
    pratio =  (Likelihoodnew/Likelihood)*ratioPrior;
end
        
alpha = min([1  pratio*proposal_ratio]); 
u = rand;                                    
if u <= alpha
    THETA = new_THETA; 
    accept = 1;
else
    accept = 0;
end
        
end

%% END
        