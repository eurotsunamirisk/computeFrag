%%% Function for Fragility estimation of Bayesian Fragility Model
%%% Written By: Hossein Ebrahimian, University of Naples Federico II (UNINA)
%%% Laters Updated: 09/2022
    
function fragility = calculate_fragility_basic(IM,theta,NDS,GRM)

if strcmp(GRM,'logit')
    Plink = @(X,bo,b1) 1./(1+exp(-bo-b1*log(X)));
elseif strcmp(GRM,'probit')
    Plink = @(X,bo,b1) normcdf(bo+b1*log(X));
elseif strcmp(GRM,'comploglog')
    Plink = @(X,bo,b1) 1-exp(-exp(bo+b1*log(X)));
end

fragility = zeros(length(IM),NDS); 
for j=1:NDS
    fragility(:,j) = Plink(IM,theta{1,j}(1),theta{1,j}(2));
end

end


