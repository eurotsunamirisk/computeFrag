%% Function for likelihood estimation of Bayesian Fragility Model
%%% written by: Hossein Ebrahimian July 2021
    
function fragility = calculate_fragility(IM,theta,NDS,GRM)

if strcmp(GRM,'logit')
    Plogit = @(X,bo,b1) (1-1./(1+exp(-bo-b1*log(X))));
elseif strcmp(GRM,'probit')
    Plogit = @(X,bo,b1) (1-normcdf(log(X/bo)/b1));
end

PDS_IM = cell(1,NDS+1); 

fragility = zeros(length(IM),NDS); 

for j=1:NDS+1
    Pl = [];
    Pl(:,1) = Plogit(IM,theta(1),theta(2));
    if j<=NDS
        for k=2:j
            Pl(:,k) = Plogit(IM,theta(2*k-1),theta(2*k)).*(1-sum(Pl,2));
        end
        PDS_IM{1,j} = Pl(:,end);
    else
        for k=2:NDS
            Pl(:,k) = Plogit(IM,theta(2*k-1),theta(2*k)).*(1-sum(Pl,2));
        end
        PDS_IM{1,j} = 1-sum(Pl,2);
    end
end

sumP = zeros(size(IM));
for j=NDS+1:-1:2
    sumP = sumP + PDS_IM{1,j};
    fragility(:,j-1) = sumP;
end

end



