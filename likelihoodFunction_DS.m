%%% Function for likelihood estimation of Bayesian Fragility Model
%%% Written By: Hossein Ebrahimian, University of Naples Federico II (UNINA)
%%% Laters Updated: 09/2022

function [y,lny] = likelihoodFunction_DS(dataDS,theta,GRM)

%% Define 1-?

if strcmp(GRM,'logit')
    Plink = @(X,bo,b1) (1-1./(1+exp(-bo-b1*log(X))));
elseif strcmp(GRM,'probit')
    Plink = @(X,bo,b1) (1-normcdf(bo+b1*log(X)));
elseif strcmp(GRM,'comploglog')
    Plink = @(X,bo,b1) (exp(-exp(bo+b1*log(X))));
end

%% Calculate P(DS|IM)

NDS = length(dataDS)-1;
PDS_IM = cell(1,NDS+1); 

for j=1:NDS+1
    Pl = [];
    Pl(:,1) = Plink(dataDS{1,j},theta(1),theta(2));
    if j<=NDS
        for k=2:j
            Pl(:,k) = Plink(dataDS{1,j},theta(2*k-1),theta(2*k)).*(1-sum(Pl,2));
        end
        PDS_IM{1,j} = Pl(:,end);
    else
        for k=2:NDS
            Pl(:,k) = Plink(dataDS{1,j},theta(2*k-1),theta(2*k)).*(1-sum(Pl,2));
        end
        PDS_IM{1,j} = 1-sum(Pl,2);
    end
end

%% Calculate the likelihood and loglikelihood

% y=1;
% for i=1:NDS+1
%     y = y*prod(PDS_IM{1,j});
% end

lny=0;
for j=1:NDS+1
    lny = lny+sum(log(PDS_IM{1,j}));
end
y = exp(lny);

end



