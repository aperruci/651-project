function [tr, r, tth, theta] = approxTrajectory(X,T,nr,nth,method,tm,rm)
%APPROXTRAJECTORY Approximate (r,theta) trajectory from generating
%function.
%   Input:
%   X      = boundary conditions [initial r (orbital radius)
%                                 initial theta (central angle)
%                                 final r
%                                 final theta                  ]
%                                 initial r_dot                ]
%                                 initial theta_dot            ]
%                                 final rdot                   ]
%                                 final theta_dot              ]
%   T      = time of flight
%   nr     = number of evaluation points for radius
%   nth    = number of evaluation points for central angle
%   method = "TH" tangent hyperbolic OR
%            "CP" cubic polynomial OR
%            "2CP" two-joined cubic polynomial (phasing only)
%   tm = segment switch time for phasing problem
%   rm = segment switch radius for phasing problem
%
%   Output:
%   tr    = time steps for radius
%   r     = approximate orbit radius over time
%   tth   = time steps for central angle
%   theta = approximate central angle over time

% Check input:
if nargin < 5
    method = "TH";
end

% Boundary conditions:
r0 = X(1);
theta0 = X(2);
r1 = X(3);
theta1 = X(4);
rdot0 = X(5);
thetadot0 = X(6);
rdot1 = X(7);
thetadot1 = X(8);

% Discretize t:
tr = linspace(0,T,nr);
tth = linspace(0,T,nth);

% Approximate trajectory:
switch method
    case "TH" % Tangent hyperbolic
        omega = 3; % shaping parameter, 1<=omega<=3
        ar = r0;
        br = r1;
        ath = theta0;
        bth = theta1;
        t0 = T/2;

        r = 0.5*((ar + br) + (br - ar)*tanh((tr - t0)/omega)); % Eqn 14
        theta = 0.5*((ath + bth) + (bth - ath)*tanh((tth - t0)/omega)); % Eqn 15
    
    case "CP" % Cubic polynomial
        a = (2*(r0 - r1) + (rdot0 + rdot1)*T)/T^3; % Eqn C3
        b = -(3*(r0 - r1) + (2*rdot0 + rdot1)*T)/T^2;
        c = rdot0;
        d = r0;
        
        e = (2*(theta0 - theta1) + (thetadot0 + thetadot1)*T)/T^3;
        f = -(3*(theta0 - theta1) + (2*thetadot0 + thetadot1)*T)/T^2;
        g = thetadot0;
        h = theta0;
        
        r = a*tr.^3 + b*tr.^2 + c*tr + d;
        theta = e*tth.^3 + f*tth.^2 + g*tth + h;
    
    case "2CP" % Two-jointed cubic polynomials
        as1 = (2*(r0 - rm) + rdot0*tm)/tm^3; % Eqn D3 is missing terms in the paper
        bs1 = -(3*(r0 - rm) + 2*rdot0*tm)/tm^2;
        cs1 = rdot0;
        ds1 = r0;

        D = (T - tm)^3;
        as2 = (2*(rm - r1) + rdot1*T - rdot1*tm)/D;
        bs2 = -(tm*(3*rm - 3*r1 + rdot1*T) + rdot1*T^2 - 2*rdot1*tm^2 - T*(3*r1 - 3*rm))/D;
        cs2 = -(rdot1*tm^3 - 2*rdot1*T^2*tm + T*tm*(6*r1 - 6*rm + rdot1*tm))/D;
        ds2 = (rm*T^3 - rdot1*T^2*tm^2 - 3*rm*T^2*tm + rdot1*T*tm^3 + 3*r1*T*tm^2 - r1*tm^3)/D;
        
        e = (2*(theta0 - theta1) + (thetadot0 + thetadot1)*T)/T^3; % Eqn C3
        f = -(3*(theta0 - theta1) + (2*thetadot0 + thetadot1)*T)/T^2;
        g = thetadot0;
        h = theta0;
        
        rs1 = as1*tr(tr<tm).^3 + bs1*tr(tr<tm).^2 + cs1*tr(tr<tm) + ds1;
        rs2 = as2*tr(tr>=tm).^3 + bs2*tr(tr>=tm).^2 + cs2*tr(tr>=tm) + ds2;
        r = [rs1, rs2];
        theta = e*tth.^3 + f*tth.^2 + g*tth + h;

    otherwise
        error("Invalid method. ''method'' must be ''TH'',''CP'', or ''2CP''")
end
        