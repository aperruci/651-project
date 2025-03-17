function [tr, r, tth, theta] = approxTrajectory(X,T,nr,nth,method,omega)
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
%   omega = width parameter for TH method. 1 <= omega <= 3
%
%   Output:
%   tr    = time steps for radius
%   r     = approximate orbit radius over time
%   tth   = time steps for central angle
%   theta = approximate central angle over time

% Check input:
if nargin < 6 && method == "TH"
    omega = 2;
end
if nargin < 5
    method = "TH";
    omega = 2;
end
if all(method ~= ["TH", "CP", "2CP"])
    error("Invalid method. ''method'' must be ''TH'',''CP'', or ''2CP''")
end
if any([omega<1,omega>3])
    error("Invalid value of omega. Must satisfy 1 <= omega <= 3")
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
        ar = r0;
        br = r1;
        ath = theta0;
        bth = theta1;
        t0 = T/2;

        r = 0.5*((ar + br) + (br - ar)*tanh((tr - t0)/omega)); % Eqn 14
        theta = 0.5*((ath + bth) + (bth - ath)*tanh((tth - t0)/omega)); % Eqn 15
    
    case "CP" % Cubic polynomial
        a = (2*(r0 - r1) + (rdot0 + rdot1)*T)/T^3;
        
end
        