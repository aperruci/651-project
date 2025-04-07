
%% LEOtoGEO.m
%
% Computes the Finite Fourier Series (FFS) approximation of constrained 
% low-thrust LEO-to-GEO rendezvous.
%
% AER E 6510 Spring 2025
% Final Project
% Alex Perruci and Hunter Underwood

clear
close all

% Parameters:
Nrev = 7;           % number of revolutions
nr = 2;             % number of radius FFS terms
ntheta = 3;         % number of theta FFS terms
Isp = 3.7183;       % specific impulse (TU)
m0 = 1e5/5.9722e24; % initial mass (MU)
nDP = 40;           % number of discretization points
Ta_max = 0.0153;    % maximum thrust acceleration (DU/TU^2)
T = 148.73;         % final time (TU)

% Boundary Conditions:
r0 = 1.0313;           % initial radius (DU)
theta0 = 0;            % initial polar angle (rad)
r1 = 6.61;             % final radius (DU)
theta1 = 47.123;       % final polar angle (rad)
rdot0 = 0;             % initial radial velocity (DU/TU)
thetadot0 = 0.9562;    % initial angular velocity (rad/TU)
rdot1 = 0;             % final radial velocity (DU/TU)
thetadot1 = 0.058842;  % final radial velocity (rad/TU)
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

% Only keep the free parameters to optimize:
X_guess = [coeffs_guess(1)            % a0
           coeffs_guess(6:2*nr+1)     % a3, b3, ..., anr, bnr
           coeffs_guess(2*nr+2)       % c0
           coeffs_guess(2*nr+7:end)]; % c3, d3, ... cntheta, dntheta

% TODO Put for loop here to iterate over weight
weight = 0;

% FFS Optimization:
t_DP = linspace(0,T,nDP);
deltaV_only = false;
[X_opt,fit] = fmincon(@(X)objectiveFun(X,t_DP,nr,ntheta,BC,Isp,m0,weight,deltaV_only), ...
                X_guess, ...
                [],[],[],[],[],[], ...
                @(X)constraintFun(X,t_DP,nr,ntheta,BC,Ta_max));
[a0_opt,a_opt,b_opt,c0_opt,c_opt,d_opt] = params2Coeffs(X_opt,t_DP,nr,ntheta,BC);

% Equation of motion residual:
f_res = objectiveFun(X_opt,t_DP,nr,ntheta,BC,Isp,m0,0,deltaV_only);

% Delta V:
deltaV_only = true;
delta_V = objectiveFun(X_opt,t_DP,nr,ntheta,BC,Isp,m0,0,deltaV_only);

% Trajectory:
t = linspace(0,T,1000);
[Ta,r,theta] = trajectoryFFS(t,a0_opt,a_opt,b_opt,c0_opt,c_opt,d_opt);
[x,y] = pol2cart(theta,r);

% TODO Plot every case
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
