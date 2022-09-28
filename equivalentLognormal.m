%%% Equivalent Lognormal Fragility Parameters
%%% Written By: Hossein Ebrahimian, University of Naples Federico II (UNINA)
%%% Laters Updated: 09/2022

function [eta,beta] = equivalentLognormal(x,fragility)

NDS = size(fragility,2);

eta = zeros(NDS,1);
x16 = zeros(NDS,1);
x84 = zeros(NDS,1);

for j=1:NDS
    eta(j) = interpola(fragility(:,j),x,0.50);
    x16(j) = interpola(fragility(:,j),x,normcdf(-1.0));
    x84(j) = interpola(fragility(:,j),x,normcdf(1.0));
end

beta = 0.5*log(x84./x16);

end