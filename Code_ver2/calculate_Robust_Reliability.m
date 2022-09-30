%%% This m-file is for Calculating Robust Reliability
%%% Written By: Hossein Ebrahimian, University of Naples Federico II (UNINA)
%%% Laters Updated: 09/2022

% fragility: size [length(IM)*(number of samples)]

%% Main
 
function [firstM,sigma] = calculate_Robust_Reliability(fragility)

firstM = mean(fragility,2);
secondM = mean(fragility.^2,2);
sigma = sqrt(secondM-firstM.^2);

end

%%