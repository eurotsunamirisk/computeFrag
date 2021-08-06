%%% This scripts uses MCMC simulation with Metropolis-Hastings procedure
%%% to estimate the probability P(theta|Data) 
%%% written by: Hossein Ebrahimian July 2021

function [samples,priorPDF_par] = posteriorFragilityFunction_logit(dataDS,dalpha0,dalpha1,THETA)

%% Initialize the Metropolis-Hastings sampler

maxIterations = 600;         % Maximum number of iterations
burnin        = 100;         % burnin period
numChain      = 5;           % number of chain

do_plot = 0;

NDS = length(dataDS)-1;

vec_threshold = 2.5;
alpha0 = -ceil(max(abs(THETA(1:2:end))))*vec_threshold:dalpha0:ceil(max(abs(THETA(1:2:end))))*vec_threshold;
alpha1 = -ceil(max(abs(THETA(2:2:end))))*vec_threshold:dalpha1:ceil(max(abs(THETA(2:2:end))))*vec_threshold;

Name_uncerParameter = {};
vec_uncerParameter = {};
for i=1:NDS
    Name_uncerParameter = [Name_uncerParameter,['\alpha_{0,',num2str(i-1),'}'],['\alpha_{1,',num2str(i-1),'}']];  % name of Uncertain Parameters
    vec_uncerParameter =  [vec_uncerParameter,alpha0,alpha1];
end

numUP = length(Name_uncerParameter);   % number of Uncertain Parameters
nbins = 40;                            % number of bins for the histogram

%% Define the PDF for Prior 

priorPDF = {};
priorPDF_par = {};
priorPDFfunction = {};

for i=1:NDS
    
    priorPDF = [priorPDF,'normal','normal'];
    priorPDF_par = [priorPDF_par,[THETA(2*i-1),0.50*abs(THETA(2*i-1))],[THETA(2*i),0.50*abs(THETA(2*i))]];
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

likelihoodFunction = @likelihoodFunction_DS_logit;

DATA = dataDS;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Start sampling

nchain = 1;    

figureName = cell(1,numChain);
    
while nchain <= numChain
    
%     display(['           - Chain Number = ',num2str(nchain)]) 
    
    state  = zeros(numUP,maxIterations);  % Storage space for our samples
    accept = zeros(numUP,maxIterations);  % Storage space for accept/reject decisions

    for iter = 1:maxIterations
        for n = 1:numUP
            if nchain==1
                [THETA,accept(n,iter)] = sample_posterior_MCMC(THETA,n,proposalPDF{1,n},proposalPDF_par{1,n},priorPDF{1,n},priorPDF_par{1,n},likelihoodFunction,DATA);
            else
                [THETA,accept(n,iter)] = sample_posterior_MCMC(THETA,n,'kernel',psample{1,n},priorPDF{1,n},priorPDF_par{1,n},likelihoodFunction,DATA);
            end
        end        
        state(:,iter) = THETA;        
    end
        
    %%% Samples we take for further analysis
    samples = state(:,burnin+1:1:maxIterations);
    
    %%% Define the Kernel density
    psample = cell(1,numUP);
    for n = 1:numUP
        psample{1,n} = fitdist((samples(n,:))','kernel');
        THETA(n,1) = psample{1,n}.random;
    end
    
    %%% Plot Samples
    if do_plot
    figure
    figureName{1,nchain} = ['Chain ',num2str(nchain),' / ',num2str(numChain)];
    set(gcf,'name',figureName{1,nchain},'numbertitle','off')
    set(gcf,'Position',[50 450-(nchain-1)*100 1750 170])
    for n = 1:numUP
    subplot(1,numUP,n)    
    h = [];
    nameh = {};
    [Xbins,PMF,h(1)] = plot_hist(samples(n,:),min(vec_uncerParameter{1,n}),max(vec_uncerParameter{1,n}),nbins,1);
    hold on
    if nchain > 1
        hold on
        [fn,xn] = ksdensity(psample{1,n}.InputData.data,'function','pdf');
        h(2) = plot(xn,fn/max(fn)*max(PMF),'-','color',[0.50,0.50,0.50],'LineWidth',1.5);
    end
    hold on
    PDFprior = priorPDFfunction{1,n}(vec_uncerParameter{1,n},priorPDF_par{1,n});
    h(3)=plot(vec_uncerParameter{1,n},PDFprior/max(PDFprior)*max(PMF)/2,'--','color',[0.9290, 0.6940, 0.1250],'LineWidth',1.50);
    hold off
    xlim([min(vec_uncerParameter{1,n}),max(vec_uncerParameter{1,n})])
    xlabel(Name_uncerParameter{1,n},'fontsize',14)
    if n==1
    ylabel('PDF','fontsize',14)
    end
    set(gca,'fontsize',12)    
    end
    end
 
    nchain = nchain + 1;
    
end

%% Plot the RV Distributions

if do_plot
for nchain=1:numChain
    close (figureName{1,nchain})
end
end

% figure
% set(gcf,'Position',[50 350 1750 170])
% for n = 1:numUP
%         
%     subplot(1,numUP,n)
%     [Xbins,~,~] = plot_hist(samples(n,:),min(vec_uncerParameter{1,n}),max(vec_uncerParameter{1,n}),nbins,'default');
%     hold on
%     priorPDF = priorPDFfunction{1,n}(Xbins,priorPDF_par{1,n});
%     plot(Xbins,priorPDF/sum(priorPDF),'-','color',[0.9290, 0.6940, 0.1250],'LineWidth',2)
%     xlim([min(vec_uncerParameter{1,n}),max(vec_uncerParameter{1,n})])
%     xlabel(Name_uncerParameter{1,n},'fontsize',20)
%     if n==1
%     ylabel('PMF','fontsize',20)
%     end
%     %legend('sample posterior by MCMC','prior')
%     set(gca,'fontsize',14)
%     
% end

%% END

