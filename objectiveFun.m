function out = objectiveFun(X,t,nr,ntheta,BC,Isp,m0)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[a0,a,b,c0,c,d] = params2Coeffs(X,t,nr,ntheta,BC);
% Equation of motion residual:
[Ta,r,~,rdot,thetadot,rddot,thetaddot] = trajectoryFFS(t,a0,a,b,c0,c,d);
mu = 1;
f = r.^2 .* (thetadot.*rddot - rdot.*thetaddot) ...
    + thetadot.*(mu - 2.*r.*rdot.^2) - (r.*thetadot).^3;
f_res = dot(f,f);

% Delta v:
g0 = 1;
g0Isp = g0*Isp;
m = zeros(length(t),1);
m(1) = m0;
mdot = zeros(length(t),1);
mdot(1) = -abs(Ta(1)*m(1) / g0Isp);
for k = 2:length(t)
    m(k) = m(k-1) + mdot(k-1)*(t(k)-t(k-1));
    mdot(k) = -abs(Ta(k-1)*m(k-1) / g0Isp);
end
deltaV = g0Isp*log(m0/m(end));

% Cost function:
out = f_res + 1e-3*deltaV;
