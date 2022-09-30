%%% This Function performs a GLR model on the data for estimating the priors
%%% Written By: Hossein Ebrahimian, University of Naples Federico II (UNINA)
%%% Laters Updated: 09/2022


function [THETA,DSj,Ncol] = calculate_glm_basic(DS,NDS,IM,GRM)

%% GLM fit

GRMCoef    = cell(1,NDS);
p_collapse = cell(1,NDS);
Ncol       = cell(1,NDS);

DSj=[];
for j=1:NDS+1
    DSj = [DSj;DS{1,j}];
end

for j=2:NDS+1
    num_collapse = [];
    for k=1:j-1
        num_collapse = [num_collapse;zeros(length(DS{1,k}),1)];
    end
    for k=j:NDS+1
        num_collapse = [num_collapse;ones(length(DS{1,k}),1)];
    end
    Ncol{1,j-1} = num_collapse;
    
    [GRMCoef{1,j-1},~] = glmfit(log(DSj),num_collapse,'binomial','link',GRM);
    p_collapse{1,j-1} = glmval(GRMCoef{1,j-1},log(IM),GRM)';
    
end

%% Find THETA

THETA = GRMCoef;

end

%% END