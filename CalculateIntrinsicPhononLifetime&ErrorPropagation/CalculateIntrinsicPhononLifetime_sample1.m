clear
clc

T = [300 250 200 80];     %temperature 
tau = [0.3332 0.3655 0.3942 0.4151];    %P3/P1
dt = [0.0039 0.0071 0.007 0.0072];       %error bar 

X = tau(1);	
dx = dt(1);
Y = tau(4);	
dy = dt(4);
T = 287.9;
avg = T/2*(log(Y)-log(X))^-1; % intrinsic phonon lifetime
% error propagation
err = sqrt( (T/2*(log(Y)-log(X))^-2*(1/X)*dx)^2 + (T/2*(log(Y)-log(X))^-2*(1/Y)*dy)^2 );
result = [avg-err avg avg+err err]