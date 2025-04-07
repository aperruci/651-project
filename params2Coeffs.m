function [a0,a,b,c0,c,d] = params2Coeffs(X,t,nr,ntheta,BC)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% Boundary conditions:
r0 = BC(1);
theta0 = BC(2);
r1 = BC(3);
theta1 = BC(4);
rdot0 = BC(5);
thetadot0 = BC(6);
rdot1 = BC(7);
thetadot1 = BC(8);
T = t(end);

% FFS Coefficients:
a0 = X(1);
a(3:nr) = X(2:2:2*nr-3);
a(1) = (r0 - r1)/2 - sum(a(3:2:end));
a(2) = (r0 + r1 - a0)/2 - sum(a(4:2:end));

b(3:nr) = X(3:2:2*nr-3);
b(1) = T/(2*pi)*(rdot0 - rdot1) - sum((3:2:nr).*b(3:2:end));
b(2) = T/(4*pi)*(rdot0 + rdot1) - 0.5*sum((4:2:nr).*b(4:2:end));

c0 = X(2*nr-3 + 1);
c(3:ntheta) = X(2*nr-3 + 2:2:end);
c(1) = (theta0 - theta1)/2 - sum(c(3:2:end));
c(2) = (theta0 + theta1 - c0)/2 - sum(c(4:2:end));

d(3:ntheta) = X(2*nr-3 + 3:2:end);
d(1) = T/(2*pi)*(thetadot0 - thetadot1) - sum((3:2:ntheta).*d(3:2:end));
d(2) = T/(4*pi)*(thetadot0 + thetadot1) - 0.5*sum((4:2:ntheta).*d(4:2:end));
end
