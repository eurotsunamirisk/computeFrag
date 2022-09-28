%%% This scripts uses MCMC simulation with Metropolis-Hastings procedure to estimate the probability P(theta|Data) 
%%% Written By: Hossein Ebrahimian, University of Naples Federico II (UNINA)
%%% Laters Updated: 09/2022

function [samples,priorPDF_par] = posteriorFragilityFunction(dataDS,dalpha0,dalpha1,THETA,MCMCtype,GRM,COVprior)

%% Initialize the Metropolis-Hastings sampler

maxIterations_chain1     = 1020;  % Maximum number of iterations
maxIterations_chain2_end = 2000;  % Maximum number of iterations
burnin        = 20;               % burnin period
numChain      = 6;                % number of chain
kernelType    = 'adaptive';       % kernelType = 'adaptive' / 'nonadaptive'

do_plot = 0;

NDS = length(dataDS)-1;

vec_threshold = 14.0; % 6 before   
% alpha0 = -ceil(max(abs(THETA(1:2:end))))*vec_threshold:dalpha0:ceil(max(abs(THETA(1:2:end))))*vec_threshold;
% alpha1 = -ceil(max(abs(THETA(2:2:end))))*vec_threshold:dalpha1:ceil(max(abs(THETA(2:2:end))))*vec_threshold;

Name_uncerParameter = {};
vec_uncerParameter = {};
for i=1:NDS
    alpha0 = -ceil(abs(THETA(2*i-1)))*vec_threshold:dalpha0:ceil(abs(THETA(2*i-1)))*vec_threshold;
    alpha1 = -ceil(abs(THETA(2*i)))*vec_threshold:dalpha1:ceil(abs(THETA(2*i)))*vec_threshold;
    Name_uncerParameter = [Name_uncerParameter,['\it\alpha\rm_{0,',num2str(i-1),'}'],['\it\alpha\rm_{1,',num2str(i-1),'}']];  % name of Uncertain Parameters
    vec_uncerParameter =  [vec_uncerParameter,alpha0,alpha1];
end

numUP = 2*NDS;      % number of Uncertain Parameters
nbins = 30;         % number of bins for the histogram

%% Define the PDF for Prior 

priorPDF = {};
priorPDF_par = {};
priorPDFfunction = {};

for i=1:NDS
    
    priorPDF = [priorPDF,'normal','normal'];
    priorPDF_par = [priorPDF_par,[THETA(2*i-1),COVprior(2*i-1)*abs(THETA(2*i-1))],[THETA(2*i),COVprior(2*i)*abs(THETA(2*i))]];
    priorPDFfunction = [priorPDFfunction,{@(x,par) normpdf(x,par(1),par(2))},{@(x,par) normpdf(x,par(1),par(2))}];
    
%     priorPDF = [priorPDF,'uniform','uniform'];
%     priorPDF_par = [priorPDF_par,{[min(alpha0),max(alpha0)],[min(alpha1),max(alpha1)]}];
%     priorPDFfunction = [priorPDFfunction,{@(x,par) unifpdf(x,par(1),par(2))},{@(x,par) unifpdf(x,par(1),par(2))}];

end

%% Proposal Distributions

proposalPDF = {};
proposalPDF_par = {};

for i=1:NDS
    proposalPDF = [proposalPDF,'normal','normal'];
    proposalPDF_par = [proposalPDF_par,0.30,0.30];    % COV
end

%% likelihoodFunction definition

likelihoodFunction = @likelihoodFunction_DS;

DATA = [dataDS,GRM];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Start sampling

nchain = 1;    

figureName = cell(1,numChain);
    
while nchain <= numChain
    
    display(['           - Chain Number = ',num2str(nchain)])
    
    if nchain==1
        maxIterations = maxIterations_chain1;
    else
        maxIterations = maxIterations_chain2_end;
    end    
    
    state  = zeros(numUP,maxIterations);  % Storage space for our samples
    accept = zeros(numUP,maxIterations);  % Storage space for accept/reject decisions

    for iter = 1:maxIterations
        display(['              --- Iteration Number = ',num2str(iter),'/',num2str(maxIterations)]) 
        
        if nchain==1
            for n = 1:numUP
                [THETA,accept(n,iter)] = sample_posterior_MCMC(THETA,n,proposalPDF{1,n},proposalPDF_par{1,n},priorPDF{1,n},priorPDF_par{1,n},likelihoodFunction,DATA);
            end
        else
            if strcmp(MCMCtype,'Blockwise')
                n=1;
                while n<=numUP
                    [THETA_dummy,accept(:,iter)] = sample_posterior_MCMC_updated(THETA,seeds,weights,priorPDF,priorPDF_par,[],likelihoodFunction,DATA);
                    if ~isequal(THETA_dummy,THETA)
                        THETA = THETA_dummy;
                        break
                    end
                    n=n+1;
                end
            elseif strcmp(MCMCtype,'Compwise')
                for n = 1:numUP
                    [THETA,accept(n,iter)] = sample_posterior_MCMC(THETA,n,'kernel',{seeds,weights},priorPDF{1,n},priorPDF_par{1,n},likelihoodFunction,DATA);
                end
            end
        end
        
        state(:,iter) = THETA;
        
    end
        
    %% Samples we take for further analysis
    
    if nchain==1
        seeds = state(:,burnin+1:maxIterations);
    else
        [~,iix] = unique(state','rows');
        seeds = state(:,sort(iix));
        %seeds = state;
    end
    
    weights = calculateWeights(seeds, kernelType);
   
end

samples = seeds;
    
end

%% END

