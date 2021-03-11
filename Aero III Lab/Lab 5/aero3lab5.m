%% Experimental Aerodynamics III Lab 5: Electric Motor Performance Analysis
% Steps
%   Load in structure and separate data
%   Calculations
%   Plotting
%% Preallocation and specifiers
% File structure: 
% | RPM (1) | Torque (N-m) (2) | Voltage (V) (3) | Current (I) (4) |
% Unclear if current units are Amps or mA or another derivation of Amp
numMots = 2; % number of tested motors
trials = 3; % number of trials per motor

%data = struct();
legSpec = {'M1 R1','M1 R2', 'M1 R3', 'M2 R1', 'M2 R2', 'M2 R3'}; % legend entries
pSpec = {'bo','g*','rv'; 'co', 'm*', 'yv'}; % plot specifier for type of data point ixj = 2x3 matrix
lSpec = {'b','g','r'; 'c','m','y'}; % plot specifier for line color for just line plot

t = ones(2,3); ip = t; op = t; me = t;
%% Loop
% Loop backwards to create entire structure array first
%   (so that it does not change size every iteration)
for i = numMots:-1:1
    for j = trials:-1:1
        file = sprintf('motor%d-run%d.dat',[i j]); % current file to load in
        data = importdata(file);
        % Temp vars
        rpm = data(:,1);
        tau = data(:,2).*0.0010; % N-mm to N-m
        volts = data(:,3);
        current = data(:,4);
        
        % Calculations
        N = rpm./60; % temp var
        
        motor(i).run(j).pout = 2.*pi.*N.*tau./1000; % kW
        motor(i).run(j).pin = volts.*current./1000; % kW
        motor(i).run(j).eta = motor(i).run(j).pout./motor(i).run(j).pin;
        
        % Allocation of rest of vars needed for plots
        motor(i).run(j).rpm = rpm;
        motor(i).run(j).tau = tau;
        
        
        % Torque vs. RPM
        figure(1)
        hold on
        t(i,j) = plot(motor(i).run(j).rpm, motor(i).run(j).tau, pSpec{i,j}, 'linewidth',1.15);
        plot(motor(i).run(j).rpm, motor(i).run(j).tau, lSpec{i,j})
        title('Torque vs. RPM')
        xlabel('RPM')
        ylabel('Torque (N-m)')
        
        % Input Power vs. RPM
        figure(2)
        hold on
        ip(i,j) = plot(motor(i).run(j).rpm, motor(i).run(j).pin, pSpec{i,j}, 'linewidth',1.15);
        plot(motor(i).run(j).rpm, motor(i).run(j).pin, lSpec{i,j})
        title('Input Power vs. RPM')
        xlabel('RPM')
        ylabel('Power In (kW)')
        
        % Shaft Output Power vs. RPM
        figure(3)
        hold on
        op(i,j) = plot(motor(i).run(j).rpm, motor(i).run(j).pout, pSpec{i,j}, 'linewidth',1.15);
        plot(motor(i).run(j).rpm, motor(i).run(j).pout, lSpec{i,j})
        title('Shaft Output Power vs. RPM')
        xlabel('RPM')
        ylabel('Power Out (kW)')
        
        % Motor Efficiency vs. RPM
        figure(4)
        hold on
        me(i,j) = plot(motor(i).run(j).rpm, motor(i).run(j).eta, pSpec{i,j}, 'linewidth',1.15);
        plot(motor(i).run(j).rpm, motor(i).run(j).eta, lSpec{i,j})
        title('Motor Efficiency vs. RPM')
        xlabel('RPM')
        ylabel('Efficiency')
    end
end
% Legend entries
phs = {t, ip, op, me}; % plot handles
% Matching the legend to each plot
for plot = 1:4 % number of plots
    ph = phs{plot}; % current plot handle
    figure(plot)
    % Format below tells legend which to plot, format based on legSpec format
    legend([ph(1,1), ph(1,2), ph(1,3), ph(2,1), ph(2,2), ph(2,3)], legSpec)
end
% Prevent coplotting for next run of code
hold off all