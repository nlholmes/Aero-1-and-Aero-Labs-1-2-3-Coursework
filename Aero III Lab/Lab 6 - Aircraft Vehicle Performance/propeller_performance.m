function outdata = propeller_performance()
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


data = importdata('Lab-3_prop_data_all_sessions_clean.dat');

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
%etap = {};
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
            diavec(i) = mean(dia(:,:));
        %end
    end
outdata.CT = C_T;
outdata.Ctau = C_tau;
outdata.CP = C_P;
outdata.etap = etap_plot;
outdata.v = vinf; % need to plot against velocity
outdata.J = J;
outdata.dia = diavec; % diameters of all props
end

