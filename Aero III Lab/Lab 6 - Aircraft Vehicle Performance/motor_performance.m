function [motor] = motor_performance()
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
            
            % Averaging the runs for aircraft performance analysis
            pout{i,j} = 2.*pi.*N.*tau./1000;
            pin{i,j} = volts.*current./1000;
            eta{i,j} = pout{i,j}./pin{i,j};
            rpmvec{i,j} = rpm;
            tauvec{i,j} = tau;
        end
        %{
        % Taking average values from each run
        motor(i).rpm = mean(rpm{i,:});
        motor(i).tau = mean(tau{i,:});
        motor(i).pout = mean(pout{i,:});
        motor(i).pin = mean(pin{i,:});
        motor(i).eta = mean(eta{i,:});
        %}
        
        % Average values
        % could try motor.rpm(i,:) = to get into mat format for further
        % calculations
        %{
        motor(i).rps = mean(cell2mat(rpmvec{i,:}))./60;
        motor(i).tau = mean(cell2mat(tauvec{i,:}));
        motor(i).pout = mean(cell2mat(pout{i,:}));
        motor(i).pin = mean(cell2mat(pin{i,:}));
        motor(i).eta = mean(cell2mat(eta{i,:}));
        %}
        motor(i).rpm = rpmvec{i,:};
        motor(i).tau = tauvec{i,:};
        motor(i).pout = pout{i,:};
        motor(i).pin = pin{i,:};
        motor(i).eta = eta{i,:};
    end
end

