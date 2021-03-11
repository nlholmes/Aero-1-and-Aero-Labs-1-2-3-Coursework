%% Experimental Aerodynamics III, Lab 3: Propeller Analysis
% NEEDS:
% Clean up
% Optional: optimize code (make loops better, try to find those unique param places 
%   like discussed below)

% Plots*******************************************************************
%   1 all curves for each prop = 7
%   (1 single curve vs J) * (4 curves) * (3 varying parameters, rest fixed) = 12 plots
%       FIX THE LOOP LIKE DID FOR PITCH
%       ADD LEGENDS, TITLES, AXES LABELS
%           Need to automate these somehow (well, thats optional for
%           robust, reusable code at least
%           Might just have to come back to this kind of thing after report
%           is written, so for now just type in legends and plot errata
%           manually
%       REGIONS OF GRAPH WITH THAT ONE MATLAB COMMAND WITH SHADING
%% Load in Data and Define Constants
% Data Format:
%   no props (1)| prop dia (in) (2)| prop pitch (in) (3)| throttle (4)| Q (psf) (5)| V (volts) (6)|
%   I (amps) (7)| Thrust (lb) (8)| Torque (in-lb) (9)| RPM (10)| Material (1-APC, 2-Wood)(11)|
% 1 in = 0.0254 m
% 1 psf = 47.88 pa
% 1 lb = 4.448 N
% 1 in-lb = 0.113 N-m


data = importdata('Lab-3_prop_data_all_sessions.dat');

% Data adjustment to include only 100% throttle, can work for other
% throttle settings too
maxThrotInds = data(:,4) == 100; % finds the cases where 100 throttle used
data = data(maxThrotInds,:); % gets portion of data only where 100 throttle used



% Constants
rho = 1.18; % freestream density, kg/m^3

%% Attempt at finding indicies and data for uniques
%{
%[uvals, ur, uc] = unique(data(:,1:3)); % the first three columns
% One for each column
[u1vals, r1] = unique(data(:,1)); % returns the unique values and the rows that they start
[u2vals, r2] = unique(data(:,2));
[u3vals, r3] = unique(data(:,3));
%}
[uvals, rind] = unique(data(:,1:3), 'rows'); % uvals are values of uniques
% 'rows' treats each row as a separate entity
% as specified data(:,1:3) it compares first 3 cols together and looks for
% unique combinations of those amongst the rows (what we want)
% rind returns the inds of the first
%   the prop for each goes from the rind int to (next rind int - 1)
% length of rind represents number of propellers tested
rind = sort(rind); % sorts rind so that indicies are in numerically increasing order

%% Finding indicies in rind to match up with their parametric study partner
% uvals tells:
%    1 to 3 for pitch
%    3 to 6 for diameter
%    3 and 7 for blade number, 7 has larger dia by 0.5 but thats the closest
%[N, edges, bin] = histcounts(uvals);
%uvals2 = unique(uvals);
%test = ismember(uvals,uvals2);

% Ideas:
%   Use unique on uvals per col
%   Find size of uvals and loop thru rows, returning row indicies where
%       one col differs, find other rows that have same two other cols
%           no need to specify that same one col, because uvals already set
%           up like that
%% Loading data for each propeller based off of unique inds
num_elems = [diff(rind)',length(data(:,1)) - rind(end)];
% diff gets the difference between adjacent, transposed b/c it rets col vec
%   this shows number of elems as rind is row inds
%   last num elem is not included as rind was instructed to find the first occurance
%   used the row length of the data matrix and the final index given by rind to find this final num elem

% Can see that the propellors do not have an equal number of data points
% Can use cells
% Nested for (i: length num_elems aka length rind, j: num_elems(i))
% Cannot plot in loop because of trade studies
%   Calculated data for plots must therefore be in cells
%       Coeffs, etap_plot, J

C_P = {};
C_T = {};
C_tau = {};
J = {};
etap = {};
etap_plot = {};
for i = 1:length(num_elems)
    %for j = num_elems(i)
        r = rind(i); % temporary rind val
        ne = num_elems(i); % temporary num_elems
        ind2 = r + ne - 1;
        % Separate data, elements per prop have (data(r:(r+n),col)
        num = data(r:ind2,1);
        dia = data(r:ind2,2) .* 0.0254; % in to m
        pitch = data(r:ind2,3) .* 0.0254; % in to m
        throt = data(r:ind2,4);
        Q = data(r:ind2,5) .* 47.88; % psf to pa
        volts = data(r:ind2,6);
        current = data(r:ind2,7);
        thrust = data(r:ind2,8) .* 4.448; % lb to N
        torque = data(r:ind2,9) .* 0.113; % in-lb to N-m
        rpm = data(r:ind2,10);
        %mat = data(r:ind2,11);
        
        
        % Calculations
        n = rpm./60; % rot/sec
        
        % Freestream velocity
        vinf = sqrt(2.*Q./rho);
        
        % Power, P = V*I
        P = volts.*current; 
        
        % Advance ratio to normalize freestream velocity
        J{i} = vinf./(n.*dia); % cell for plot
        
        % Coefficients, in cells for plot
        C_T{i} = thrust./(rho.*n.^2.*dia.^4);
        C_tau{i} = torque./(rho.*n.^2.*dia.^5);
        C_P{i} = P./(rho.*n.^3.*dia.^5);
        
        % Propeller efficiency
        etap = C_T{i}.*J{i}./(2.*pi.*C_tau{i}); 
        % Setting etap to zero where torque or thrust is negative
        etap_fix = etap; % new variable representing the to be corrected etap
        etap_fix(torque < 0 | thrust < 0) = 0; % correction for etap
        etap_plot{i} = etap_fix; % cell for plot
        
    %end
end


%% Plots for each Prop
% Plot
        %   Trade Studies (4 per study = 4*3 = 12)
        %       Pitch, diameter, blade number as constants
        %   Efficiencies (1 per prop = 7)
        %   Total number of plots = 19
% So plotting multiple throttle settings gives multiple 'curves'
%   Need to adjust data to only include 100% throttle...?

% Need to label regions for these plots
propReg = cell(1,7);
abrakeReg = propReg;
millReg = abrakeReg;
for t = 1:length(rind)
    %thisProp = 7; % propeller to plot, nums 1 thru 7
    figure(t) % coplot all for each propeller
    % Use patches to find regimes
    hold on
    plot(J{t}, etap_plot{t},'rv')
    plot(J{t}, C_T{t}, 'bo')
    plot(J{t}, C_tau{t}.*10, 'k.','markersize',10) % multiplied by 10
    plot(J{t}, C_P{t}, 'r*')
    
    % Using patches for regions
    %   Propeller region (GREEN): thrust and torque > 0
    propReg{t} = C_T{t} > 0 & C_tau{t} > 0;
    propRange = J{t}(propReg{t}); % J values of where prop region condition satisfied
    propxr = [min(propRange) max(propRange) max(propRange) min(propRange)];
    propyr = [-0.4 -0.4 1 1];
    patch(propxr,propyr,'green','FaceAlpha',0.2);
    
    %   Air brake region (CYAN): thrust < 0 and torque > 0
    abrakeReg{t} = C_T{t} < 0 & C_tau{t} > 0;
    abrakeRange = J{t}(abrakeReg{t}); % J values of where prop region condition satisfied
    % abrakexr should have minimim value of max(propRange) to eliminate white space
    abrakexr = [max(propRange) max(abrakeRange) max(abrakeRange) max(propRange)];
    abrakeyr = [-0.4 -0.4 1 1];
    patch(abrakexr,abrakeyr,'cyan','FaceAlpha',0.2);
    
    %   Windmill region (BLACK): thrust and torque < 0
    millReg{t} = C_T{t} < 0 & C_tau{t} < 0;
    millRange = J{t}(millReg{t}); % J values of where prop region condition satisfied
    % millxr should have min val of max(abrakeRange) to elim white space
    millxr = [max(abrakeRange) max(millRange) max(millRange) max(abrakeRange)];
    millyr = [-0.4 -0.4 1 1];
    patch(millxr,millyr,'black','FaceAlpha',0.2);
    
    
    legend('\eta_p','C_T','C_\tau (x10)','C_P')
    xlabel('J')
    ylabel('\eta_p, C_T, C_\tau, C_P')
    title(sprintf('Plots for Propeller %d',t))
    hold off
end



%% Parametric studies, fig nums dep on last val of t
%    1 to 3 for pitch (all 4 curves coplotted)
%    3 to 6 for diameter (all 4 curves coplotted)
%    3 and 7 for blade number, 7 has larger dia by 0.5 but thats the
%       closest (all 4 curves coplotted)

% Below are the vectors representing which propellers have the indicated varying parameter
pitch = [1 2 3]; % param 1
blade_dia = [3 4 5 6]; % param 2
blade_num = [3 7]; % param 3
% Could nest another loop in here and use ifs to choose which one to plot
% to make a little more efficient/neater I think
%for p = 1:maxLengthParams % change this it is ruining plots
%    figcount = t+1; % start for figure count
num_params = 3;
figcount = t + 1; % starting figure count
for p = 1:num_params    
    % Pitch (1)
    % Blade Diameter (2)
    % Blade Number (3)
    if p == 1 % if the current param to plot is the first param (pitch)
        figcount = plotPropParams(J, etap_plot, C_P, C_T, C_tau, pitch, figcount);
    elseif p == 2 % blade_dia
        figcount = plotPropParams(J, etap_plot, C_P, C_T, C_tau, blade_dia, figcount);
    else % p == 3, plot blade_num
        figcount = plotPropParams(J, etap_plot, C_P, C_T, C_tau, blade_num, figcount);
    end
end

%{
%% Old separated data
start = 1; % orig value 16
stop = 32; % orig value 32

num = data(start:stop,1);
dia = data(start:stop,2) .* 0.0254; % in to m
pitch = data(start:stop,3) .* 0.0254; % in to m
throt = data(start:stop,4);
Q = data(start:stop,5) .* 47.88; % psf to pa
volts = data(start:stop,6);
current = data(start:stop,7);
thrust = data(start:stop,8) .* 4.448; % lb to N
torque = data(start:stop,9) .* 0.113; % in-lb to N-m
rpm = data(start:stop,10);
mat = data(start:stop,11);

%% Old Calculations

% Rotational speed
n = rpm./60; % rot/sec
             
% Freestream velocity             
vinf = sqrt(2.*Q./rho);
% Power, P = V*I
P = volts.*current; 
% Advance ratio to normalize freestream velocity
J = vinf./(n.*dia); 
% Coefficients
C_T = thrust./(rho.*n.^2.*dia.^4);
C_tau = torque./(rho.*n.^2.*dia.^5);
C_P = P./(rho.*n.^3.*dia.^5);
% Propeller efficiency
etap = C_T.*J./(2.*pi.*C_tau); 
% Setting etap to zero where torque or thrust is negative
etap_fix = etap; % new variable representing the to be corrected etap
etap_fix(torque < 0 | thrust < 0) = 0; % correction for etap
%}