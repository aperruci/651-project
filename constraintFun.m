function [C,Ceq] = constraintFun(X,t,nr,ntheta,BC,Ta_max)
%CONSTRAINTFUN Summary of this function goes here
%   Detailed explanation goes here

[a0,a,b,c0,c,d] = params2Coeffs(X,t,nr,ntheta,BC);
Ta = trajectoryFFS(t,a0,a,b,c0,c,d);
C = (Ta./Ta_max).^2 - 1;
Ceq = 0; % trivial constraint for fmincon
