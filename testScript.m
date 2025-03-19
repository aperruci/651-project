%% testScript
% Use to debug functions.

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

% Initial trajectory:
[tr_cp, r_cp, tth_cp, theta_cp] = approxTrajectory(X,T,nr,nth,'CP');
[tr_th, r_th, tth_th, theta_th] = approxTrajectory(X,T,nr,nth,'TH');
figure
tiledlayout(2,1,'tilespacing','tight','padding','tight')
nexttile
plot(tr_cp,r_cp,...
     tr_th,r_th, ...
     'LineWidth',2)
title('Initial trajectory')
legend('cubic polynomial','tanh')
ylabel('r (DU)')
nexttile
plot(tth_cp,theta_cp, ...
     tth_th,theta_th, ...
     'LineWidth',2)
xlabel('t (TU)')
ylabel('\theta (rad)')