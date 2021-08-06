%%% This Function performs a GLR model on the data for estimating the priors
%%% written by: Hossein Ebrahimian


function THETA = calculate_glm(DS,NDS,IM,GRM)

%% Simple logistic regression to see the data (Lj)

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
    
    if strcmp(GRM,'logit')
        
        [GRMCoef{1,j},~] = glmfit(log(DSj),num_collapse,'binomial','link',GRM);
        p_collapse{1,j} = glmval(GRMCoef{1,j},log(IM),GRM)';
    
    elseif strcmp(GRM,'probit')
        
        [GRMCoef_probit,~] = glmfit(log(DSj),[num_collapse,ones(length(num_collapse),1)],'binomial','link',GRM);
        b0 = GRMCoef_probit(1);
        b1 = GRMCoef_probit(2);
        
        % convert probit coefficients to lognormal distribution parameters
        GRMCoef{1,j} = [exp(-b0/b1);1/b1];
        p_collapse{1,j} = normcdf(log(IM),log(GRMCoef{1,j}(1)),GRMCoef{1,j}(2));
        
    end
end

%% Find THETA

THETA = reshape(cell2mat(GRMCoef),2*NDS,1);

end

%% END