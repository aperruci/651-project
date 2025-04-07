function out = objectiveFun(X,t,nr,ntheta,BC,Isp,m0,weight,deltaV_only)
% OBJECTIVEFUN Objective function for FFS coefficient optimization.
% Enforces the equations of motion and penalizes delta V.
%
% Input:  X      = unknown FFS coefficients
%         t      = time vector
%         nr     = number of radius FFS terms
%         ntheta = number of theta FFS terms
%         BC     = boundary conditions
%         Ta_max = maximum thrust acceleration (DU/TU^2)
%         weight = multiplier for delta V in the fitness
%
% Output: out    = fitness

% Equation of motion residual:
[a0,a,b,c0,c,d] = params2Coeffs(X,t,nr,ntheta,BC);
[Ta,r,~,rdot,thetadot,rddot,thetaddot] = trajectoryFFS(t,a0,a,b,c0,c,d);
mu = 1; % standard gravitational parameter (DU^3/TU^2)
f = r.^2 .* (thetadot.*rddot - rdot.*thetaddot) ... % equation of motion 
    + thetadot.*(mu - 2.*r.*rdot.^2) - (r.*thetadot).^3;
f_res = dot(f,f); % equation of motion should be 0

% Delta v:
g0 = 1; % gravity acceleration (DU/TU^2)
g0Isp = g0*Isp;
m = zeros(length(t),1); % mass (MU)
m(1) = m0;
mdot = zeros(length(t),1); % mass flow rate (MU/TU)
mdot(1) = -abs(Ta(1)*m(1) / g0Isp);
for k = 2:length(t) % iteratively calculate mass
    m(k) = m(k-1) + mdot(k-1)*(t(k)-t(k-1));
    mdot(k) = -abs(Ta(k-1)*m(k-1) / g0Isp);
end
deltaV = g0Isp*log(m0/m(end)); % (DU/TU)

% Fitness:
if deltaV_only
    out = deltaV;
else
    out = f_res + weight*deltaV;
end