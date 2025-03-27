%% LEO-GEO Maneuver
clear
close all

% Parameters:
Nrev = 7;
nr = 2;
ntheta = 3;
Isp = 3.7183;
m0 = 1e5/1.9891e30;
nDP = 40;
Ta_max = 0.0153;
T = 148.73;

% Boundary Conditions:
r0 = 1.0313;
theta0 = 0;
r1 = 6.61;
theta1 = 47.123;
rdot0 = 0;
thetadot0 = 0.9562;
rdot1 = 0;
thetadot1 = 0.058842;
BC = [r0
      theta0
      r1   
      theta1
      rdot0
      thetadot0
      rdot1
      thetadot1];

% Initial guess for FFS coefficients:
coeffs_guess = initialCoeffs(BC,T,nr,ntheta,'CP');
X_guess = [coeffs_guess(1)            %a0
           coeffs_guess(6:2*nr+1)     %a3-anr, b3_nr
           coeffs_guess(2*nr+2)       %c0
           coeffs_guess(2*nr+7:end)]; %c3-cntheta, d3-dntheta

% FFS Optimization:
t_DP = linspace(0,T,nDP);
[X_opt,fit] = fmincon(@(X)objectiveFun(X,t_DP,nr,ntheta,BC,Isp,m0), ...
                X_guess, ...
                [],[],[],[],[],[], ...
                @(X)constraintFun(X,t_DP,nr,ntheta,BC,Ta_max));

a0_opt = X_opt(1);
a_opt(3:nr) = X_opt(2:2:2*nr-3);
a_opt(1) = (r0 - r1)/2 - sum(a_opt(3:2:end));
a_opt(2) = (r0 + r1 - a0_opt)/2 - sum(a_opt(4:2:end));

b_opt(3:nr) = X_opt(3:2:2*nr-3);
b_opt(1) = T/(2*pi)*(rdot0 - rdot1) - sum((3:2:nr).*b_opt(3:2:end));
b_opt(2) = T/(4*pi)*(rdot0 + rdot1) - 0.5*sum((4:2:nr).*b_opt(4:2:end));

c0_opt = X_opt(2*nr-3 + 1);
c_opt(3:ntheta) = X_opt(2*nr-3 + 2:2:end);
c_opt(1) = (theta0 - theta1)/2 - sum(c_opt(3:2:end));
c_opt(2) = (theta0 + theta1 - c0_opt)/2 - sum(c_opt(4:2:end));

d_opt(3:ntheta) = X_opt(2*nr-3 + 3:2:end);
d_opt(1) = T/(2*pi)*(thetadot0 - thetadot1) - sum((3:2:ntheta).*d_opt(3:2:end));
d_opt(2) = T/(4*pi)*(thetadot0 + thetadot1) - 0.5*sum((4:2:ntheta).*d_opt(4:2:end));

% Trajectory:
t = linspace(0,T,1000);
r = 0.5*a0_opt + sum(a_opt'.*cos((1:nr)'*pi/T.*t) + b_opt'.*sin((1:nr)'*pi/T.*t));
rdot = sum(-a_opt' .* ((1:nr)'*pi/T) .* sin((1:nr)'*pi/T.*t) ...
           + b_opt' .* ((1:nr)'*pi/T) .* cos((1:nr)'*pi/T.*t));
rddot = sum(-a_opt' .* ((1:nr)'*pi/T).^2 .* cos((1:nr)'*pi/T.*t) ...
            - b_opt' .* ((1:nr)'*pi/T).^2 .* sin((1:nr)'*pi/T.*t));      
theta = 0.5*c0_opt + sum(c_opt'.*cos((1:ntheta)'*pi/T.*t) + d_opt'.*sin((1:ntheta)'*pi/T.*t));
thetadot = sum(-c_opt' .* ((1:ntheta)'*pi/T) .* sin((1:ntheta)'*pi/T.*t) ...
               + d_opt' .* ((1:ntheta)'*pi/T) .* cos((1:ntheta)'*pi/T.*t));
thetaddot = sum(-c_opt' .* ((1:ntheta)'*pi/T).^2 .* cos((1:ntheta)'*pi/T.*t) ...
                - d_opt' .* ((1:ntheta)'*pi/T).^2 .* sin((1:ntheta)'*pi/T.*t));
[x,y] = pol2cart(theta,r);

% Thrust acceleration:
cosalpha = r.*thetadot./sqrt(rdot.^2 + (r.*thetadot).^2);
Ta = (2*rdot.*thetadot + r.*thetaddot) ./ cosalpha;

% Plots
figure
plot(x,y,'k')
xlabel('X (DU)')
ylabel('Y (DU)')

figure
plot(t,Ta,'k')
yline(Ta_max,'--k')
yline(-Ta_max,'--k')
xlabel('t (TU)')
ylabel('T_a (DU/TU^2)')