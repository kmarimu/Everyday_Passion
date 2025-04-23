function [x3,p3] = simulatePolicy3(x0,z,pChemDrive,a,tau,etac,etad,pcMax,xMax,xMin,t,hDeadline,xStar)
% simulatePolicy3 simulates the third electric vehicle charging policy
% (when plugged in and below minimum charge, charge just fast enough to
% meet a deadline).
%
% Inputs:
%   x0, the battery's initial chemical energy, kWh
%   z, a K x 1 of indicators that the vehicle is plugged in
%   pChemDrive, a K x 1 of chemical powers discharged to drive, kW
%   a, a scalar discrete-time dynamics parameter
%   tau, the battery's self-dissipation time constant, h
%   etac, the battery's charging efficiency
%   etad, the battery's discharging efficiency
%   pcMax, the battery's maximum charging electrical power, kW
%   xMax, the battery's chemical energy capacity, kWh
%   xMin, the minimum acceptable chemical energy, kWh
%   t, the simulation time span, h
%   hDeadline, the hour of day of the charging deadline (0 = midnight)
%   xStar, the desired charge at the deadline, kWh
%
% Outputs:
%   x3, a K+1 x 1 vector of stored chemical energies, kWh
%   p3, a K x 1 vector of electrical charging powers, kW

% timing
K = length(z); % number of time steps
dt = t(2) - t(1); % time step duration, h

% initialization
x3 = zeros(K+1,1); % stored chemical energy, kWh
x3(1) = x0; % initial state
pChem3 = -pChemDrive; % chemical charging power, kW (initialized with power used for driving)
y3 = zeros(K,1); % indicator of charging mode

% simulation
for k=1:K
    % charging decision
    if z(k) == 0 || x3(k) == xMax
        y3(k) = 0;
    elseif z(k) == 1 && x3(k) < xMin
        y3(k) = 1;
    end
    if y3(k) == 1
        kStar = mod(hDeadline - mod(t(k), 24), 24);
        stepsRemaining = max(1, (kStar / dt));
        sum_a_terms = (1 - a^stepsRemaining) / (1 - a);
        pChem3(k) = min(etac * pcMax, (xStar - a^stepsRemaining * x3(k)) / ((1 - a) * sum_a_terms * tau));
        y3(k+1) = y3(k);
    end
    
    % dynamic update (applies to all cases)
    x3(k+1) = a*x3(k) + (1-a)*tau*pChem3(k);
end
p3 = max(pChem3/etac,etad*pChem3); % charging electrical power, kW
p3(z==0) = 0; % no electrical charging or discharging while unplugged

end

