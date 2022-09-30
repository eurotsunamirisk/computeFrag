function [rfragility, sfragility, etaIMc, betaIMc,...
    sample_theta_model] =...
    RF_y_logReg_model(damage_data, NDS, D, IM,...
    dvec_alpha0, dvec_alpha1, confidence)
% K. Trevlopoulos
% August 2021
% This is a super-function, which is used to avoid MATLAB cells and
% structures as input

GRM = 'logit';

DS = cell(1,NDS+1);
for j=1:NDS+1
    DS{1,j} = damage_data( damage_data(:,2) == D(j), 1);
end

theta_prior_model = calculate_glm(DS, NDS, IM, GRM);

[sample_theta_model, priorPDF_par_model] =...
    posteriorFragilityFunction_logit(...
    DS, dvec_alpha0, dvec_alpha1, theta_prior_model);

RF_model = RobustFragility_DS(...
    sample_theta_model, IM, NDS, confidence, GRM);

rfragility = RF_model.fragility;
sfragility = RF_model.sfragility;
etaIMc = RF_model.etaIMpl;
betaIMc = RF_model.betaIMpl;

end

