function [outdata] = engine_performance()
%% Constants
fuelrho = 862.75; % fuel density, kg/m^3
r = 0.4445; % length of lever, m
Qhv = 11.3.*1e6; % heating value of fuel, MJ/kg to J/kg
etac = 0.99; % combustion efficiency
g = 9.81; % accel due to gravity, m/s^2, used for force calculation


%% Load in data
% File format
% RPM 1 (1) | RPM 2 (2) | RPM 3 (3) | RPM 4 (4) | RPM 5 (5)| 
% Force 1 (gms) (6) | Force 2 (gms) (7) | Force 3 (gms) (8) | Force 4 (gms) (9) |Force 5 (gms) (10)|
% Force Correction (gms) (11)| Time for 10cc fuel consumption (seconds) (12)

% Use structure array this time (for each lab session)
data = importdata('Lab-4_Tuesday-Session-Data_20191015.dat');
%data = importdata('Lab-4_Wednesday-Session-1-Data_20191016.dat');
rpm = data(:,1:5)./2; % divide by 2 as the data is for 2 prop blades
force_raw = data(:,6:10).*1e-3; % gms (grams) to kgm
force_corr = data(:,11).*1e-3; % gms to kgm
fuel_time = data(:,12); % number of seconds for 10ccs to go down
force = (force_raw - force_corr).*g; % mult by gravity to get real force, Newtons (kgm to N)

% SHOULD THIS BE MULT? (same then as line below it)
fuel = 10e-6./fuel_time; % fuel consumption, 10cc / number of seconds it took, 1e6 cc = 1 m^3, units m^3/s
%fuel = fuel_time.*10e-6;
%fuel = 10e-6./fuel_time.^(-1);
%   there are 10e-6 cubic meters in 10 cc
%{
%% Structure array for data
% filenames is a cell array variable containing names of all files
filenames = {'Lab-4_Tuesday-Session-Data_20191015.dat',...
    'Lab-4_Wednesday-Session-1-Data_20191016.dat',...
    'Lab-4_Wednesday-Session-2-Data_20191016.dat',...
    'Lab-4_Thursday-Session-Data_20191017.dat'};
for i = 1:length(filenames) % for number of files to load in
    filedata = importdata(filenames{i}); % current file to load in
    
    data(i).rpm = filedata(:,1:5)./2; % divide by 2 as the data is for 2 prop blades
    data(i).force_raw = filedata(:,6:10).*1e-3; % gms (grams) to kg
    data(i).force_corr = filedata(:,11).*1e-3; % gms to kg
    data(i).fuel_time = filedata(:,12); 
    data(i).force = (data(i).force_raw - data(i).force_corr).*g; % mult by gravity to get real force, Newtons
    
    fuel = 10e-6./data(i).fuel_time;
end
%}
%% Calculations
% Averaging data for calculations
avgrpm = mean(rpm,2); % mean by rows
avgforce = mean(force,2);

% Torque
tau = avgforce.*r;

% Power
    % RPS of shaft
    n = avgrpm./60;
P = 2.*pi.*n.*tau; % watts

% Thermal efficiency
%   Mass Flow Rate = density * volumetric flow rate
    mdot = fuelrho.*fuel; % kg/s
    %mdot = fuelrho./fuel;
etath = P./(mdot.*Qhv.*etac);
%etath = 1e3.*P./(mdot.*Qhv.*etac); % using kW instead of watts but now has dimensions

outdata.tau = tau;
outdata.P = P;
outdata.rps = n;
outdata.etath = etath;
end

