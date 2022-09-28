%%%   Script for calculating the Robust Fragility for empirical fragility estimation and its confidence band 
%%% Written By: Hossein Ebrahimian, University of Naples Federico II (UNINA)
%%% Laters Updated: 09/2022

function RF = RobustFragility_DS(sample_theta,x,NDS,confidence,GRM,tol)

%% Sampling Procedure

sample_fragility = cell(1,1);
sample_PDS_IM    = cell(1,1);
dummy = [];
dstep = 5;

count=0;
for i=1:size(sample_theta,2)
    [dummy_fragility,dummy_PDS] = calculate_fragility(x',sample_theta(:,i),NDS,GRM);
    if ( any(any(diff(dummy_fragility(1:dstep:length(x),:))<tol)>0) || any(dummy_fragility(1,:)>0.01) ) 
        dummy = [dummy,i];
    else
        count=count+1;
        sample_fragility{1,count} = dummy_fragility;
        sample_PDS_IM{1,count} = cell2mat(dummy_PDS);
    end
end

rejected_samples = dummy;
sample_theta(:,rejected_samples) = [];

%% Calculate the Robust Fragility and its standard deviation

Rfragility = zeros(length(x),NDS);
sfragility = zeros(length(x),NDS);
sample_fragility_DSj = cell(1,NDS);

for j=1:NDS
    sample_fragility_DSj{1,j} = zeros(length(x),size(sample_theta,2));
    for i=1:size(sample_theta,2)
        sample_fragility_DSj{1,j}(:,i) = sample_fragility{1,i}(:,j);
    end
    [Rfragility(:,j),sfragility(:,j)] = calculate_Robust_Reliability(sample_fragility_DSj{1,j});
end

%% Save

RF = struct;

RF.fragility            = Rfragility;
RF.sfragility           = sfragility;
RF.sample_fragility     = sample_fragility;
RF.sample_fragility_DSj = sample_fragility_DSj;
RF.rejected_samples     = rejected_samples;
RF.sample_theta         = sample_theta;

RF.RFp = Rfragility+confidence*sfragility;
RF.RFm = Rfragility-confidence*sfragility;

RF.RF84 = Rfragility+sfragility;
RF.RF16 = Rfragility-sfragility;

for j=1:NDS
    RF.etaIMc(j)= interpola(Rfragility(:,j),x,0.50);
    RF.IMc16(j) = interpola(Rfragility(:,j),x,normcdf(-1.0));
    RF.IMc84(j) = interpola(Rfragility(:,j),x,normcdf(1.0));
    
    RF.etaIMc_RF84(j) = interpola(RF.RF84(:,j),x,0.50);
    RF.etaIMc_RF16(j) = interpola(RF.RF16(:,j),x,0.50);
end

RF.betaIMc = 0.5*log(RF.IMc84./RF.IMc16);
RF.betaUF  = 0.5*log(RF.etaIMc_RF16./RF.etaIMc_RF84);

RF.sample_PDS_IM = sample_PDS_IM;

end

%% END