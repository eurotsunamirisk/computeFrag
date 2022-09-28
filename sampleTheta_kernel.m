%%% Function for sampling from a Kernel
%%% Written By: Hossein Ebrahimian, University of Naples Federico II (UNINA)
%%% Laters Updated: 09/2022

%% Main

function sample_theta = sampleTheta_kernel(seeds,weights)

F = cumsum(weights./sum(weights));
iKernel = find(F>=rand,1,'first');

Mu = seeds(:,iKernel);
S = (weights(iKernel))^2*cov(seeds');
L = chol(S)';             % L = lower triangular matrix L*L'=S. Matlab chol function gives the upper triangular matrix L'
sample_theta = Mu+L*randn(size(Mu));
 
return

%% END