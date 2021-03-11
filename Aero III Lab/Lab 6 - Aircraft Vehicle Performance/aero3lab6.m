%% Experimental Aerodynamics 6: Aircraft Performance Analysis - the final lab of college
% Make a function for this
% Inputs:
%   Operating velocities                    V
%   Zero-lift drag coefficient              CD0
%   Wing area                               S
%   Take-off mass                           m1
%   Oswalled Efficiency Factor              e
%   Aspect Ratio                            AR
%   Flight time                             t_flight
%   Electric motor or gas engine?           'Motor' or 'Engine' or a derivative of either or 'gas' or 'elec'
% Additional Inputs -- Data:
%   Propeller
%   Electric motor/engine

% Constants:
%   Density
%   Fuel heating value
%   Acceleration due to gravity
%   Combustion Efficiency

% Outputs:
%   Try to figure out plot handles and figure handles and see if can feed those through
%   Best E, R, and ROC
%% Calculations
q = 1/2.*rho.*V.^2.; % Dynamic Pressure
PR = q.*V.*CD0.*S + W1.^2.*V./(q.*S).*(1./(pi.*e.*AR)); % Power Required
PA = etap.*Pavs; % Power Available, need to find this using data, process detailed below
% Pavs (P_available_source) is either P_eng or P_mot
% P is shaft power available
% To calculate it, first assume engine is operating at peak power
% RPM at peak power used for operating range advanced ratios
J = V./(n.*D); % n is RPS for peak power and D is prop diameter
% Use a 3rd degree polynomial fit for the eta_prop (etap) vs. J curve to
% evaluate the etap across the entire flight range
% Use polyfit

% Slope of the minimum drag line
beta = arctand(min(PR./V)); % ?, slope is arctan of (PR/V)|min
% Tangential point of min drag line (MDL) and PR curve gives:
%   Maximum theoretical range
% Minimum of PR curve gives (point of minimum PR):
%   Maximum theoretical endurance

% Breguet equations for power rated aircrafts at a constant angle-of-attack 
% and constant altitude are used to evaluate the endurance and range
E = etap./Cp.*sqrt(2.*rho.*S).*CL.^(3/2)./CD.*(1./sqrt(W2) - 1./sqrt(W1)); % Endurance
% Cp is power specific fuel consumption

% W2 is final weight of aircraft
Cp = mDotf.*g./Pavs; % mdotf =  mass flow rate of fuel for the peak engine/motor power
% When it can not be measured, it can be equivalently calculated:
mDotf = Pin./(Qhv.*etac); % Pin is input power
% FOR THE CASE OF AN ELECTRIC MOTOR
% mDotf is not real fuel, so it has an equaivalent
% Pin_elec_equiv = mDotf*Qhv*etac
%   do not have etac nor Qhv vals nor mDotf vals, this just shows that if
%   compared to the engine equation for mDotf, that Pin_elec_equiv = mDotf
Pin_elec_equiv = V*I; % real variable name will just be Pin
%   Voltage times Current

% Lift and drag coefficients
CL = W1./(q.*S);
CD = CD0 + CL.^2./(pi.*e.*AR);

% Final aircraft weight
W2 = W1 - t_flight.*mDotf.*g; % t_flight is total flight time

% Range for all operating velocities
R = etap./Cp.*CL./CD.*log(W1./W2);
ROC = (PA - PR)./W1; % rate of climb

% Best E, R, and ROC are the maximum values of these curves
%   will want to find the velocity at which they occur
%% Calling Functions for previous lab codes
g = 9.81;

airframe = 1; % can also = 2
type = 'gas';
propNum = 4; % 1 to 7 I think



% have nested loops, one for eng/mot type, one for numProps
numProps = 7;
type = {'gas','elec1','elec2'};
numEM = length(type);
figNum = 1; % starting figure number

numAirframes = 2; % total number of airframes to test
% For each airframe
for airframe = 1:numAirframes
    % type is engine type
    % have really gas, elec1, elec2
    % prop number is the number prop to test with that motor
    
    % Data for each airframe
    if airframe == 1
        V = [1:40];
        CD0 = 0.01;
        S = 0.5;
        m1 = 3;
        e = 0.8;
        AR = 3;
        t = 15*60;
    elseif airframe == 2
        V = [1:60];
        CD0 = 0.018;
        S = 0.8;
        m1 = 10;
        e = 0.6;
        AR = 6;
        t = 30*60;
    end
    
    % For each engine/motor
    for j = 1:numEM
        % For each propeller
        thisType = type{j};
        for i = 1:numProps
            propNum = i;
            % i is prop number below in aircraft_performance function
            [E, R, ROC, PA, PR, V, n ,CL, CD] = aircraft_performance(V,CD0,S,m1,e,AR,t,thisType,propNum);

            % Theoretical values
            % Calculating Beta and max theo range and max theo endurance locations
            minPRV = min(PR./V); % theoretical
            VminPRV = V(PR./V == minPRV); % location of max theo range, % max theoretical range (tangential point with PR curve)
            minPR = min(PR); 
            VminPR = V(PR == min(PR)); % max theoretical endurance location (location of min PR)
            
            %aminPRV = min(PR)./V(PR == min(PR)); % actual
            %aVminPRVloc = V(min(PR)./V(PR == min(PR)) == min(PR)./V);
            
            beta = atan(minPRV); % will be the same for each airframe as PR no change with prop inputs
            % line: y = m*x + b, so find b
            b = minPR - beta*(VminPR);
            minDline = beta.*V + b; % slope*x + b
            % doesnt exactly hit likely due to inaccuracies in computation
            % could also be due to experimental innacuracies that the curve does not
            % increase fast enough. The minDline does pass through the minimum however.

            % Calculating best end, range, ROC (experimental best values)
            bestE(i,j) = max(E);
            bestR(i,j) = max(R);
            bestROC(i,j) = max(ROC);
            bEloc(i,j) = V(E == bestE(i,j));
            bRloc(i,j) = V(R == bestR(i,j));
            bROCloc(i,j) = V(ROC == bestROC(i,j));
            
            
            figure(figNum)
            hold on
            plot(V,PR)
            plot(V,PA)
            plot(V,minDline)
            xline(VminPRV,'k--','linewidth',2); % max theo range loc
            xline(VminPR,'k','linewidth',2); % max theo end loc
            xline(bEloc(propNum,j),'r-.'); % best end
            xline(bRloc(propNum,j),'b-.'); % best range
            xline(bROCloc(propNum,j),'m-.'); % best ROC
            
            % Set title for each figure
            titleStr = sprintf('Airframe %d, Eng/Motor %d, Propeller %d',airframe,j,i);
            title(titleStr)
            
            legend('Power Required','Power Available','Min Drag Line',...
                'Max Theo Range','Max Theo End','Max End','Max Range','Max ROC');
            hold off
            
            figNum = figNum + 1; % updating figure number
        end
        
        %{
        % Now find which eng/mot/prop has the most bests of the best (use sum)
        %   Could make this more rigorous by ensuring the values are within
        %   acceptable ranges
        locBBE = find(bestE == max(bestE)); % velocities for bests of best
        %locBBR = V(bestR == max(bestR));
        %locBBR = max(bRloc);
        locBBR = find(bestR == max(bestR));
        locBBROC = find(bestROC == max(bestROC));
        
        % TODO *********************************************************
        % Find where bestE( experimental ) - theoretical bestE is min
        % to get theoretical bestE, use intersection points 
        
        
        % best propellers for engine in categories E, R, ROC in col vec
        bestEng = [locBBE ; locBBR ; locBBROC]; 
        %}
    end
end

% Finding minimum difference between theo and exp E,R and max ROC for each
% As velocity is what is changing, the differences in velocities (locs)
% were found 
minDiffE = min(min(abs(VminPR - bEloc)));
minDiffR = min(min(abs(VminPRV - bRloc)));
% and bestROC = bestROC
bestMaxROC = max(ROC);

mdeLoc = bEloc(abs(bEloc-VminPR) == minDiffE);
mdeLoc = mdeLoc(1);
mdrLoc = bRloc(abs(bRloc-VminPRV) == minDiffR);
mdrLoc = mdrLoc(1);
bmrocLoc = bROCloc(bestROC == bestMaxROC);


bestCombo = [mdeLoc;mdrLoc;bmrocLoc]; % ultimate best E;R;ROC velocities/locations



%{
%% Plotting
figure(1)
hold on
plot(V,PR)
plot(V,PA)
plot(V,minDline)
xline(VminPRV,'k--','linewidth',2); % max theo range loc
xline(VminPR,'k','linewidth',2); % max theo end loc
xline(bEloc(propNum),'r-.'); % best end
xline(bRloc(propNum),'b-.'); % best range
xline(bROCloc(propNum),'m-.'); % best ROC
legend('Power Required','Power Available','Min Drag Line',...
    'Max Theo Range','Max Theo End','Max End','Max Range','Max ROC');
hold off
%}
%% NEEDS
% electric stuff, correcting etap (does for both atm and idk if should even)
% both: had Pin = Pavs in aircraft_performance, might need to be PA, check
% gc and handout/lab slides

% NEW:
% Might need to fix the electric motor stuff, or it could just not be
%   powerful enough for the given airframes (it seems to work for a couple
%   on the first airframe, just poor performance, as expected as they dont
%   give as much thrust I don't think)

% Need to find the best performing combo per airframe, might be done by
%   visual inspection if necessary

% Need to plot these dudes, probably just pick the best one to plot and a
%   few others to explain reasoning

% Need to find way to distinguish between the two elec motors, probably use
%   string contain and the numbers 1 and 2, use antoher if statement to
%   change the variable indexing into the motor_performance outputs

% Need to construct loop to run through all prop and engne/motor combos to
%   see what gives best endurance, roc, range, etc
%       might not be effective means of judging if the numbers are super
%       high but out of range

% Want CL vs. CD plots of all? nah that shit the same for each airfram
