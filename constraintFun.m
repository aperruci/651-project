function [C,Ceq] = constraintFun(X,t,nr,ntheta,BC,Ta_max)
%CONSTRAINTFUN Summary of this function goes here
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

% Radius and derivatives:
r = 0.5*a0 + sum(a' .* cos((1:nr)'*pi/T.*t) + b' .* sin((1:nr)'*pi/T.*t));
rdot = sum(-a' .* ((1:nr)'*pi/T) .* sin((1:nr)'*pi/T.*t) ...
           + b' .* ((1:nr)'*pi/T) .* cos((1:nr)'*pi/T.*t));

% Theta derivatives:
thetadot = sum(-c' .* ((1:ntheta)'*pi/T) .* sin((1:ntheta)'*pi/T.*t) ...
               + d' .* ((1:ntheta)'*pi/T) .* cos((1:ntheta)'*pi/T.*t));
thetaddot = sum(-c' .* ((1:ntheta)'*pi/T).^2 .* cos((1:ntheta)'*pi/T.*t) ...
                - d' .* ((1:ntheta)'*pi/T).^2 .* sin((1:ntheta)'*pi/T.*t));

% Thrust acceleration:
cosalpha = r.*thetadot./sqrt(rdot.^2 + (r.*thetadot).^2);
Ta = (2*rdot.*thetadot + r.*thetaddot) ./ cosalpha;
C = (Ta./Ta_max).^2 - 1;
Ceq = 0; % trivial constraint for fmincon
