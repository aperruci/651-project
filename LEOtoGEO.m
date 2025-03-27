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
[a0_opt,a_opt,b_opt,c0_opt,c_opt,d_opt] = params2Coeffs(X_opt,t_DP,nr,ntheta,BC);

% Trajectory:
t = linspace(0,T,1000);
[Ta,r,theta] = trajectoryFFS(t,a0_opt,a_opt,b_opt,c0_opt,c_opt,d_opt);
[x,y] = pol2cart(theta,r);

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
