%%% This is the function for claculating the evidence
%%% Written By: Hossein Ebrahimian, University of Naples Federico II (UNINA)
%%% Laters Updated: 09/2022

function [log_evidence,meanlogLikelihood,meanlogratioP ] = calculate_logE(sample_theta,likelihoodFunction,data,priorPDF_par,priorPDF_type,GRM)

logLikelihood = zeros(1,length(sample_theta));
posteriorPDF = zeros(1,length(sample_theta));
priorPDF = zeros(1,length(sample_theta));

%% Calculate Loglikelihood and Posterior

weights = calculateWeights(sample_theta,'adaptive');
prior_par = cell2mat(priorPDF_par');

for j=1:size(sample_theta,2)
    [~,logLikelihood(j)] = likelihoodFunction(data,sample_theta(:,j),GRM);
    posteriorPDF(j) = calculateKernel(sample_theta(:,j),sample_theta,weights);
    priorPDF(j) = calculateMVN(sample_theta(:,j),prior_par(:,1),prior_par(:,2),eye(length(sample_theta(:,j))),priorPDF_type);
end

%% Calculating the prior

ratioP = posteriorPDF./priorPDF;

%% Estimate the log_evidence

meanlogLikelihood = mean(logLikelihood);

meanlogratioP = mean(log(ratioP));

log_evidence = meanlogLikelihood-meanlogratioP;

end

%% END


