function [Ta,r,theta,rdot,thetadot,rddot,thetaddot] = trajectoryFFS(t,a0,a,b,c0,c,d)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Parameters:
nr = length(a);
ntheta = length(c);
T = t(end);

% Radius and derivatives:
r = 0.5*a0 + sum(a'.*cos((1:nr)'*pi/T.*t) + b'.*sin((1:nr)'*pi/T.*t));
rdot = sum(-a' .* ((1:nr)'*pi/T) .* sin((1:nr)'*pi/T.*t) ...
           + b' .* ((1:nr)'*pi/T) .* cos((1:nr)'*pi/T.*t));
rddot = sum(-a' .* ((1:nr)'*pi/T).^2 .* cos((1:nr)'*pi/T.*t) ...
            - b' .* ((1:nr)'*pi/T).^2 .* sin((1:nr)'*pi/T.*t));  

% Central angle and derivatives:
theta = 0.5*c0 + sum(c'.*cos((1:ntheta)'*pi/T.*t) + d'.*sin((1:ntheta)'*pi/T.*t));
thetadot = sum(-c' .* ((1:ntheta)'*pi/T) .* sin((1:ntheta)'*pi/T.*t) ...
               + d' .* ((1:ntheta)'*pi/T) .* cos((1:ntheta)'*pi/T.*t));
thetaddot = sum(-c' .* ((1:ntheta)'*pi/T).^2 .* cos((1:ntheta)'*pi/T.*t) ...
                - d' .* ((1:ntheta)'*pi/T).^2 .* sin((1:ntheta)'*pi/T.*t));

% Thrust acceleration:
cosalpha = r.*thetadot./sqrt(rdot.^2 + (r.*thetadot).^2);
Ta = (2*rdot.*thetadot + r.*thetaddot) ./ cosalpha;
end
