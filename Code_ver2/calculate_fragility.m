%%% Function for fragility estimation of Bayesian Fragility Model
%%% Written By: Hossein Ebrahimian, University of Naples Federico II (UNINA)
%%% Laters Updated: 09/2022
    
function [fragility,PDS_IM] = calculate_fragility(IM,theta,NDS,GRM)

%% Define 1-pi

if strcmp(GRM,'logit')
    Plink = @(X,bo,b1) (1-1./(1+exp(-bo-b1*log(X))));
elseif strcmp(GRM,'probit')
    Plink = @(X,bo,b1) (1-normcdf(bo+b1*log(X)));
elseif strcmp(GRM,'comploglog')
    Plink = @(X,bo,b1) (exp(-exp(bo+b1*log(X))));
end

%% Calculate P(DS|IM)

PDS_IM = cell(1,NDS+1); 

fragility = zeros(length(IM),NDS); 

for j=1:NDS+1
    Pl = [];
    Pl(:,1) = Plink(IM,theta(1),theta(2));
    if j<=NDS
        for k=2:j
            Pl(:,k) = Plink(IM,theta(2*k-1),theta(2*k)).*(1-sum(Pl,2));
        end
        PDS_IM{1,j} = Pl(:,end);
    else
        for k=2:NDS
            Pl(:,k) = Plink(IM,theta(2*k-1),theta(2*k)).*(1-sum(Pl,2));
        end
        PDS_IM{1,j} = 1-sum(Pl,2);
    end
end

%% Calculate the fragility

sumP = zeros(size(IM));
for j=NDS+1:-1:2
    sumP = sumP + PDS_IM{1,j};
    fragility(:,j-1) = sumP;
end

end



