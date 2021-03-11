function [E, R, ROC, PA, PR, V, n, CL, CD] = aircraft_performance(V,CD0,S,m1,e,AR,t,type,propNum)
% Inputs
%   v           range of velocities
%   CD0         zero-lift drag coefficient
%   S           wing reference area
%   m1          takeoff mass
%   e           oswald prop efficiency
%   AR          aspect ratio
%   t           flight time
%   type        'gas' engine or 'elec' motor
%   edata       data from the engine/motor experiment
%   pdata       propeller data
% Outputs
%   E           max endurance
%   R           max range
%   ROC         max rate of climb

propData = propeller_performance();
etap = propData.etap{propNum};
J = propData.J{propNum}; % n is RPS for peak power and D is prop diameter
D = propData.dia(propNum); % prop diameter is the same for all runs of htat propeller
% Propulsive generator
if(strcmp('gas',type) == 1) % do engine analysis
    % gas engine function for those calculations
    genData = engine_performance();
    % Propeller
    %V = propData.v; % Velocity

    Pavs = max(genData.P); % assume engine is operating at peak power
    n = genData.rps(genData.P == Pavs); % rps corresponding to peak power
    %{
    if length(etap) < length(Pavs)
        Pavs = Pavs(1:length(etap));
        n = n(1:length(etap));
        V = V(1:length(etap));
    elseif length(Pavs) < length(etap)
        etap = etap(1:length(Pavs));
        J = J(1:length(Pavs));
    end
    %}
elseif(contains(type,'elec') == 1) % do electric motor analysis
    % Determine which motor number it is based off of input
    if strcmp(type,'elec1') == 1
        motNum = 1; % which motor number to test
    elseif strcmp(type,'elec2') == 1
        motNum = 2;
    end
    
    % electric motor function for those calculations
    genData = motor_performance();
    %genData.pin
    %Pavs = mean(genData(1).pin{propNum,:});
    Pavs = max(genData(motNum).pin).*1000; % this ones in kW
    %rps = genData(1).rpm{propNum,:}./60;
    rps = genData(motNum).rpm./60;
    n = rps(genData(motNum).pin.*1000 == Pavs); % rps corresponding to peak power
end
%% Creating polyfit for etap vs. J to have full range of velocities tested
% polynomial fit of 3rd degree etap vs. J
% then calculate adavance ratios for new range of velocities
%   uses rps at peak power (n in this function)
% then calculate new etap's from the polyfit by plugging in the new advance ratios
pjfit = polyfit(J,etap, 3); % polyfit for Power Available
newJ = V./(n.*D);
newEtap = polyval(pjfit,newJ); % new etap values for entire velocity range

%% Uses newEtap from here on out
% Constants
g = 9.81;
rho = 1.18;
etac = 0.99;
Qhv = 11.3e6; % J/kg fuel heating value
% Doing gas engine first
W1 = m1*g;
% Calculate E, R, ROC, necessary data for curves
q = 1/2.*rho.*V.^2.; % Dynamic Pressure

PR = q.*V.*CD0.*S + W1.^2.*V./(q.*S).*(1./(pi.*e.*AR)); % Power Required
%PR = 1/2.*rho*V.^3.*CD0.*S + W1.^2./(1/2.*rho.*V.^2.*S).*(1./(pi.*e.*AR));
% newEtap used
% may need TO ADD CORRECTION FOR LOWER THAN ZERO AND HIGHER THAN 0.99 or 1
% for ELEC MOTORS
PA = newEtap.*Pavs; % Power Available, need to find this using data, process detailed below
badInd = etap < 0 | etap > 1; % these two lines maybe fix but dunno
%PA = PA(~badInd); % values of where good PAs are
etap(badInd) = 0;

% Use a 3rd degree polynomial fit for the eta_prop (etap) vs. J curve to
% evaluate the etap across the entire flight range
% Use polyfit

% Slope of the minimum drag line
% When it can not be measured, it can be equivalently calculated:
%if type == 'gas'
%    Pin = Pavs; % think for gas this is true
%elseif type == 'elec'
    % FOR THE CASE OF AN ELECTRIC MOTOR
    % mDotf is not real fuel, so it has an equaivalent
    % Pin_elec_equiv = mDotf*Qhv*etac
    %   do not have etac nor Qhv vals nor mDotf vals, this just shows that if
    %   compared to the engine equation for mDotf, that Pin_elec_equiv = mDotf
%    Pin= V*I; % real variable name will just be Pin
    %   Voltage times Current
%end
Pin = Pavs; % I DO NOT KNOW WHICH ONE THIS IS SUPPOSED TO BE
%Pin = PA;
mDotf = Pin./(Qhv.*etac); % Pin is input power
Cp = mDotf.*g./Pavs; % mdotf =  mass flow rate of fuel for the peak engine/motor power
%beta = arctand(min(PR./V)); % ?, slope is arctan of (PR/V)|min

% Lift and drag coefficients
CL = W1./(q.*S);
CD = CD0 + CL.^2./(pi.*e.*AR);
W2 = W1 - t.*mDotf.*g; % t_flight is total flight time




% Use newEtap values for entire velocity range
E = newEtap./Cp.*sqrt(2.*rho.*S).*CL.^(3/2)./CD.*(1./sqrt(W2) - 1./sqrt(W1)); % Endurance
R = newEtap./Cp.*CL./CD.*log(W1./W2); % Range for all operating velocities
ROC = (PA - PR)./W1; % rate of climb
end

