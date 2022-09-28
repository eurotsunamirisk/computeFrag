%%% Function for Calculating ordered statistics
%%% Written By: Hossein Ebrahimian, University of Naples Federico II (UNINA)
%%% Laters Updated: 09/2022

%% Main

function [xmean,x50,x16,x84,x02,x98] = ordered_statistic(x)

xmean = mean(x,2);

x50 = median(x,2);

N = size(x,2);

x02 = (interp1((1:N)',sort(x,2)',normcdf(-2)*N,'nearest'))';
x16 = (interp1((1:N)',sort(x,2)',normcdf(-1)*N,'nearest'))';
x84 = (interp1((1:N)',sort(x,2)',normcdf(1)*N,'nearest'))';
x98 = (interp1((1:N)',sort(x,2)',normcdf(2)*N,'nearest'))';   

end

%% END




