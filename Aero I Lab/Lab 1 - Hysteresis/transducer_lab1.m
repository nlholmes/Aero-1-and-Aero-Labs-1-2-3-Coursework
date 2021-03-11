% Loading data and sorting in ascending order
incdat = sort(load('inc.dat'),'ascend'); 
decdat = sort(load('dec.dat'),'ascend');

% Extracting manometer heights and adjusting for initial height reading and
% adjusting to correct units (inches to meters, 1 inch = .0254 m)
hinc = 0.0254*(incdat(:,2) - incdat(1,2)); % data already sorted so that zero reading is first element
hdec = 0.0254*(decdat(:,2) - decdat(1,2));
% Calculating Manometer Pressure (rho*g*h, where rho is the density of water in SI and g is the gravitational constant)
rho = 997.71; % kg/m^3
g = 9.81; % m/s^2
pcalc = @(h) rho.*g.*h; % anonymous function for below calculations
ipress = pcalc(hinc);
dpress = pcalc(hdec);

% Extracting transducer pressure readings and converting from psf to N/m^2
%   1psf = 47.88026 Pa
transi = 47.88026*(incdat(:,1));
transd = 47.88026*(decdat(:,1));

%% Pressure Plots
% Plotting manometer pressure vs. sensor pressure
figure(1)
plot(ipress,transi,'ro')
hold on
plot(dpress,transd,'bo')
% Lines of best fit
linfiti = polyfit(ipress,transi,1);
linfitd = polyfit(dpress,transd,1);
fitLinei = polyval(linfiti,ipress); % polyval for applying it to the data
fitLined = polyval(linfitd,dpress);
% Plotting those lines and graph details
plot(ipress,fitLinei,'r--')
plot(dpress,fitLined,'b--')
legend('Increasing manometer data','Decreasing manometer data','Increasing line of best fit',...
    'Decreasing line of best fit');
title('Manometer pressure vs. Transducer pressure')
ylabel('Transducer Pressure (Pa)');
xlabel('Manometer Pressure (Pa)');

% Plotting difference between increasing and decreasing as compared to
%   transducer readings and plotting (mano - trans)
figure(2) % new figure for differences
hold on
diffi = abs(ipress - transi); % difference in increasing (manometer - transducer), absolute value for difference
diffd = abs(dpress - transd); % difererence in decreasing
plot(ipress,diffi,'r');
plot(dpress,diffd,'b');
title('Difference in Pressure vs. Manometer pressure');
xlabel('Manometer Pressure (Pa)');
ylabel('Pressure Difference (Pa)');
diffifit = polyfit(ipress,diffi,1);
diffdfit = polyfit(dpress,diffd,1);
diffiLine = polyval(diffifit,ipress);
diffdLine = polyval(diffdfit,dpress);
plot(ipress,diffiLine,'r--');
plot(dpress,diffdLine,'b--');
legend('Increasing Difference','Decreasing Difference','Inc. Diff. Fit Line','Dec. Diff. Fit Line');
%% Plotting difference between increasing and decreasing pressures from manometer alone
% shows human error
diffmano = abs(ipress - dpress(1:20)); % only first 20 elements of dec. pressure values to match the 20 inc. vals
figure(3)
plot(ipress,diffmano);
xlabel('Manometer pressure (Pa)');
ylabel('Difference in manometer pressure (Pa)');
title('Difference in manometer pressure vs. Manometer Pressure');
legend('Difference in Inc. and Dec.');
%% Sensor Current vs. Manometer Pressure
% Getting values for current from data file, with current in mA
% (1e-3 mA = 1 A)
incAmps = (incdat(:,3)); % mA in the third column
decAmps = (decdat(:,3));

% Current vs. Pressure plot
figure(4)
hold on
plot(ipress,incAmps,'ro');              % Increasing data plot
plot(dpress,decAmps,'bo');              % Decreasing data plot
% Fit lines and plotting
ampfiti = polyfit(ipress,incAmps,1);    % Increasing fit
ampfitd = polyfit(dpress,decAmps,1);    % Decreasing fit
ampLinei = polyval(ampfiti,ipress);     % Increasing line
ampLined = polyval(ampfitd,dpress);     % Decreasing line
plot(ipress,ampLinei,'r--');
plot(dpress,ampLined,'b--');
% Graph details
ylabel('Sensor Current (mA)');
xlabel('Manometer Pressure (Pa)');
title('Sensor Current vs. Manometer Pressure');
legend('Increasing Data','Decreasing Data','Increasing fit','Decreasing fit')
%% Extra Credit: Sensor Current vs. Transducer Pressure Readings
% Transducer pressure vs sensor current plot
figure(5)
plot(transi,incAmps,'r*')
hold on
plot(transd,decAmps,'b*')
% Fit lines and plotting
transampfiti = polyfit(transi,incAmps,1);
transampfitd = polyfit(transd,decAmps,1);
plot(transi,polyval(transampfiti,transi),'r--');
plot(transd,polyval(transampfitd,transd),'b--');
% Transampfiti and transampfitd have the same slope, which means that the
% transducer is correctly calibrated. There is only slight deviation
% between the lines of best fit, that's why the blue dotted line seems to
% almost completely cover the red dotted line.

% Graph details
xlabel('Transducer Pressure (Pa)');
ylabel('Sensor Current (mA)')
title('Sensor Current vs. Transducer Pressure');
legend('Increasing data','Decreasing data','Increasing fit line','Decreasing fit line');




