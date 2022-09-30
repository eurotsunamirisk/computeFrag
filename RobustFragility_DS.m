%   ----------------------------------------------------------------------
%   Script for calculating the Robust Fragility and its confidence band 
%   -----------------------------------------------------------------------
%   Date: July 2021
%   Writen by: Hossein Ebrahimian

function RF = RobustFragility_DS(sample_theta,x,NDS,confidence,GRM)

%% Sampling Procedure

sample_fragility = cell(1,size(sample_theta,2));

for i=1:size(sample_theta,2)
    sample_fragility{1,i} = calculate_fragility(x',sample_theta(:,i),NDS,GRM);
end

%% Calculate the Robust Fragility and its standard deviation

Rfragility = zeros(length(x),NDS);
sfragility = zeros(length(x),NDS);
sample_fragility_DSj = cell(1,NDS);

for j=1:NDS
    sample_fragility_DSj{1,j} = zeros(length(x),size(sample_theta,2));
    for i=1:size(sample_theta,2)
        sample_fragility_DSj{1,j}(:,i) = sample_fragility{1,i}(:,j);
    end
    Rfragility(:,j) = mean(sample_fragility_DSj{1,j},2);
    sfragility(:,j) = std(sample_fragility_DSj{1,j},1,2);
end

%% Save

RF = struct;

RF.fragility            = Rfragility;
RF.sfragility           = sfragility;
RF.sample_fragility     = sample_fragility;
RF.sample_fragility_DSj = sample_fragility_DSj;

RF.RFp84 = Rfragility+confidence*sfragility;
RF.RFp16 = Rfragility-confidence*sfragility;

for j=1:NDS
    RF.etaIMpl(j) = interpola(Rfragility(:,j),x,0.50);
    RF.IMp16(j)   = interpola(Rfragility(:,j),x,normcdf(-1.0));
    RF.IMp84(j)   = interpola(Rfragility(:,j),x,normcdf(1.0));
    
    RF.etaIM_RFp84(j) = interpola(RF.RFp84(:,j),x,0.50);
    RF.etaIM_RFp16(j) = interpola(RF.RFp16(:,j),x,0.50);
end

RF.betaIMpl = 0.5*log(RF.IMp84./RF.IMp16);
RF.betaUF   = 0.5*log(RF.etaIM_RFp16./RF.etaIM_RFp84);

end

%% END