%% testScript
% Use to debug functions.
clear
close all

% Test [tr, r, tth, theta] = approxTrajectory(X,T,nr,nth,method,omega,tm,rm,rdotm)
% Earth-Mars Transfer
% Boundary conditions: (all units in DU and TU)
r0 = 1;
theta0 = 0;
r1 = 1.5234;
theta1 = 9.831;
rdot0 = 0;
thetadot0 = 1;
rdot1 = 0;
thetadot1 = 0.5318;
X = [r0
     theta0
     r1
     theta1
     rdot0
     thetadot0
     rdot1
     thetadot1];

% Parameters:
Nrev = 1;
nr = 2;
nth = 5;
T = 13.447;

% Initial trajectory and FFS coefficients:
[r_coeffs_cp,th_coeffs_cp,tr_cp, r_cp, tth_cp, theta_cp] = initialCoeffs(X,T,nr,nth,'CP');
[r_coeffs_th,th_coeffs_th,tr_th, r_th, tth_th, theta_th] = initialCoeffs(X,T,nr,nth,'TH');

% Initial FFS trajectory:
t = linspace(0,T);
r_ffs_cp = 0.5*r_coeffs_cp(1) + sum(r_coeffs_cp(2:2:end).*cos((1:nr)'*pi/T*t) + r_coeffs_cp(3:2:end).*sin((1:nr)'*pi/T*t));
r_ffs_th = 0.5*r_coeffs_th(1) + sum(r_coeffs_th(2:2:end).*cos((1:nr)'*pi/T*t) + r_coeffs_th(3:2:end).*sin((1:nr)'*pi/T*t));
theta_ffs_cp = 0.5*th_coeffs_cp(1) + sum(th_coeffs_cp(2:2:end).*cos((1:nth)'*pi/T*t) + th_coeffs_cp(3:2:end).*sin((1:nth)'*pi/T*t));
theta_ffs_th = 0.5*th_coeffs_th(1) + sum(th_coeffs_th(2:2:end).*cos((1:nth)'*pi/T*t) + th_coeffs_th(3:2:end).*sin((1:nth)'*pi/T*t));

figure
tiledlayout(2,1,'tilespacing','tight','padding','tight')
nexttile
plot(tr_cp,r_cp,...
     t,r_ffs_cp,...
     tr_th,r_th, ...
     t,r_ffs_th, ...
     'LineWidth',2)
title('Earth-Mars initial trajectory')
legend('cubic polynomial','FFS CP','tanh','FFS TH')
ylabel('r (DU)')
nexttile
plot(tth_cp,theta_cp, ...
     t,theta_ffs_cp, ...
     tth_th,theta_th, ...
     t, theta_ffs_th, ...
     'LineWidth',2)
xlabel('t (TU)')
ylabel('\theta (rad)')

% LEO 90 deg Phasing
% Boundary conditions: (all units in DU and TU)
r0 = 1.0313;
theta0 = 0;
r1 = 1.0313;
theta1 = 7.8539;
rdot0 = 0;
thetadot0 = 0.9548;
rdot1 = 0;
thetadot1 = 0.9548;
X = [r0
     theta0
     r1
     theta1
     rdot0
     thetadot0
     rdot1
     thetadot1];

% Switching conditions:
tm = 6;
rm = 1.05;

% Parameters:
Nrev = 1;
nr = 3;
nth = 6;
T = 8.924;

% Initial Trajectory:
[r_coeffs_2cp,th_coeffs_2cp,tr_2cp, r_2cp, tth_2cp, theta_2cp] = initialCoeffs(X,T,nr,nth,'2CP',tm,rm);
t = linspace(0,T);
r_ffs_2cp = 0.5*r_coeffs_2cp(1) + sum(r_coeffs_2cp(2:2:end).*cos((1:nr)'*pi/T*t) + r_coeffs_2cp(3:2:end).*sin((1:nr)'*pi/T*t));
theta_ffs_2cp = 0.5*th_coeffs_2cp(1) + sum(th_coeffs_2cp(2:2:end).*cos((1:nth)'*pi/T*t) + th_coeffs_2cp(3:2:end).*sin((1:nth)'*pi/T*t));

figure
tiledlayout(2,1,'tilespacing','tight','padding','tight')
nexttile
plot(tr_2cp,r_2cp,...
     t,r_ffs_2cp, ...
     'LineWidth',2)
title('LEO Phasing initial trajectory')
legend('two-jointed cubic polynomial','FFS 2CP')
ylabel('r (DU)')
nexttile
plot(tth_2cp,theta_2cp, ...
    t,theta_ffs_2cp, ...
     'LineWidth',2)
xlabel('t (TU)')
ylabel('\theta (rad)')
