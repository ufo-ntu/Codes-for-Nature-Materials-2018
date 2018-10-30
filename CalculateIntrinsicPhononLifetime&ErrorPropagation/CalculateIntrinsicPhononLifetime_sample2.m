clear
clc

T = [300 250 80];    %temperature 
tau = [0.4028 0.4224 0.4733];   %P3/P1
dt = [0.0083 0.0243 0.0237];     %error bar 

X = tau(1);	
dx = dt(1);
Y = tau(3);	
dy = dt(3);
T = 141.4;
avg = T/2*(log(Y)-log(X))^-1;% intrinsic phonon lifetime
% error propagation
err = sqrt( (T/2*(log(Y)-log(X))^-2*(1/X)*dx)^2 + (T/2*(log(Y)-log(X))^-2*(1/Y)*dy)^2 );
result = [avg-err avg avg+err err]
