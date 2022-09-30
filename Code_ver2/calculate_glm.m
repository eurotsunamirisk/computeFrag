%%% This Function performs a GLR model on the data for estimating the priors
%%% Written By: Hossein Ebrahimian, University of Naples Federico II (UNINA)
%%% Laters Updated: 09/2022

function THETA = calculate_glm(DS,NDS,IM,GRM)

%% GLM fit

GRMCoef    = cell(1,NDS);
p_collapse = cell(1,NDS);

for j=1:NDS
    
    DSj = DS{1,j};
    num_collapse = zeros(length(DS{1,j}),1);
    
    for k=j+1:NDS+1
        DSj = [DSj;DS{1,k}];
        num_collapse = [num_collapse;ones(length(DS{1,k}),1)];
    end
    
    DSN{1,j} = DSj;
    Ncol{1,j} = num_collapse;
    
    [GRMCoef{1,j},~] = glmfit(log(DSj),num_collapse,'binomial','link',GRM);
    p_collapse{1,j} = glmval(GRMCoef{1,j},log(IM),GRM)';
    
end

%% Find THETA

THETA = reshape(cell2mat(GRMCoef),2*NDS,1);

end

%% END