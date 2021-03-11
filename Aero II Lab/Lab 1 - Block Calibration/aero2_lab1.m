%% Experimental Aerodynamics II: Supersonic Wind Tunnel Block Calibration

%% Load in data
calibs = {'400','600','800','1000','1200','1400', ...
    '1600','1800','2000','2200','2400','2600'};

% Constants and variables
gamma = 1.4;
Patm = 1:length(calibs);
P01 = Patm;
P02 = Patm;
P1 = Patm;
isentropM = Patm;
rayleighM= Patm;

dataSet = {};
% Column definitions:
% time (s) | P_01 (psi) | P_1 (psi) | P_02 (psi) | Patm (psi) | T_01 (degree F)
for i = 1:length(calibs)
    % Load data
    filename = sprintf('calibration_%s.txt',calibs{i}); % string concatenation
    dataSet{i} = importdata(filename);
    thisData = dataSet{1,i}.data;
    % Getting only values of time >= 2.5s
    timeFixDex = thisData(:,1) >= 2.5; % logical array of what rows to extract
    thisData = thisData(timeFixDex,:);
    
    % Do calculations
    % Pull pressures and convert from psi to Pa
    Patm(i) = mean(thisData(:,5))*6894.76;
    P01(i) = mean(thisData(:,2))*6894.76 + Patm(i);
    P1(i) = mean(thisData(:,3))*6894.76 + Patm(i);
    P02(i) = mean(thisData(:,4))*6894.76 + Patm(i);
    
    
    syms M
    isentrop = (1 + ((gamma-1)/2)*M^2)^(gamma/(gamma-1)) == P01(i)/P1(i);
    %isentrop = (1 + (1.2)*M^2)^(7/2) == 0.1037;
    rayleigh = ((gamma+1)^2*M^2/(4*gamma*M^2-2*(gamma-1)))^(gamma/(gamma-1))...
    *(1-gamma+2*gamma*M^2/(gamma+1)) == P02(i)/P1(i);

    isentropM(i) = double(vpasolve(isentrop,M,[0,5])); % guess between 0 and 5 to eliminate negative answers
    rayleighM(i) = double(vpasolve(rayleigh,M,[0,5]));
end

%% Plot
blockNumsRaw = 400:200:2600;
isentropMRaw = isentropM;
rayleighMRaw= rayleighM;


blockNums = blockNumsRaw(3:end);
isentropM = isentropM(3:end);
rayleighM = rayleighM(3:end);

% Ignoring first 2 data points

% Best fit lines
numtest = 4; % maximum degree polynomial to test
isenFit = [];
rayFit = isenFit;
xVals = linspace(min(blockNums),max(blockNums),100);
error = {};
percDiff = error;
meanError = error;
meanPercDiff = error;
h = [];
for k = 2:numtest
    % Polynomial fitting
    isenFit{k} = polyfit(blockNums, isentropM, k);
    rayFit{k} = polyfit(blockNums, rayleighM, k);
    isenVals = polyval(isenFit{k}, xVals);
    rayVals = polyval(rayFit{k}, xVals);
    error{k} = abs(isenVals - rayVals);
    
    h{k - 1} = figure(k - 1);
    hold on
    grid on
    plot(blockNumsRaw, isentropMRaw, 'bo','linewidth',1.25,'markersize',5); % isentropic data
    plot(blockNumsRaw, rayleighMRaw, 'ro','linewidth',1.25,'markersize',5); % rayleigh data
    errorbar(xVals, isenVals, error{k},'capsize',0);
    errorbar(xVals, rayVals, error{k},'capsize',0);
    plot(xVals,isenVals,'b') % isentropic fit
    plot(xVals,rayVals,'r') % rayleigh fit
    legend('Isentropic Data','Rayleigh Data');
    xlabel('Block Number');
    ylabel('M_\infty');
    
    meanError{k} = mean(error{k});
    percDiff{k} = error{k}./(0.5*isenVals(k)+rayVals(k)).*100;
    meanPercDiff{k} = mean(percDiff{k});
end
grid on
xlabel('Block Number');
ylabel('M_\infty');

% Polynomials on plots
figure(1)
 % power of 2
icap1 = sprintf('y = %.2s*x^2 + %.2s*x + %.2s',isenFit{2}(1),isenFit{2}(2), isenFit{2}(3));
rcap1 = sprintf('y = %.2s*x^2 + %.2s*x + %.2s',rayFit{2}(1),rayFit{2}(2),rayFit{2}(3));
text(50,1.35,icap1,'color','b')
text(50,1.15,rcap1,'color','r')
title('Quadratic Fit')

figure(2)
 % power of 3
icap2 = sprintf('y = %.2s*x^3 + %.2s*x^2 +%.2s*x + %.2s',isenFit{3}(1),isenFit{3}(2),isenFit{3}(3),isenFit{3}(4));
rcap2 = sprintf('y = %.2s*x^3 + %.2s*x^2 +%.2s*x + %.2s',rayFit{3}(1),rayFit{3}(2),rayFit{3}(3),rayFit{3}(4));
text(50,1.35,icap2,'color','b')
text(50,1.15,rcap2,'color','r')
title('Cubic Fit')

figure(3)
 % power of 4
icap3 = sprintf('y = %.2s*x^4 + %.2s*x^3 + %.2s*x^2 +%.2s*x + %.2s',isenFit{4}(1),isenFit{4}(2),isenFit{4}(3),isenFit{4}(4),isenFit{4}(5));
rcap3 = sprintf('y = %.2s*x^4 + %.2s*x^3 + %.2s*x^2 +%.2s*x + %.2s',rayFit{4}(1),rayFit{4}(2),rayFit{4}(3),rayFit{4}(4),rayFit{4}(5));
text(50,1.35,icap3,'color','b')
text(50,1.15,rcap3,'color','r')
title('Quartic Fit')
