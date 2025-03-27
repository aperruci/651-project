function [C,Ceq] = constraintFun(X,t,nr,ntheta,BC,Ta_max)
% CONSTRAINTFUN Checks the maximum thrust acceleration constraint of the
% trajectory defined by the FFS coefficients and boundary conditions.
%
% Input:  X      = unknown FFS coefficients
%         t      = time vector
%         nr     = number of radius FFS terms
%         ntheta = number of theta FFS terms
%         BC     = boundary conditions
%         Ta_max = maximum thrust acceleration (DU/TU^2)
%
% Output: C   = maximum thrust acceleration inequality constraint
%         Ceq = trivial equality constraint

[a0,a,b,c0,c,d] = params2Coeffs(X,t,nr,ntheta,BC); % all FFS coefficients
Ta = trajectoryFFS(t,a0,a,b,c0,c,d); % thrust acceleration (DU/TU^2)
C = (Ta./Ta_max).^2 - 1; % maximum thrust constraint
Ceq = 0; % trivial constraint for fmincon
