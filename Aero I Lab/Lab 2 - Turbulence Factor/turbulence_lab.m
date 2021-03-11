dat1 = load('tuesrun1.dat'); % tues run 1
dat2 = load('tuesrun2.dat'); % tues run 2
dat3 = load('wed1run1.dat'); % wed1 run 1
dat4 = load('wed1run2.dat'); % wed1 run 2
dat5 = [dat1(1,:);load('wed2run1.dat')]; % wed2 run 1, concatenated because is missing first element
dat6 = [dat1(1,:);load('wed2run2.dat')]; % wed2 run 2
dat7 = load('thurrun1.dat'); % thur run 1
dat8 = [dat1(1,:);load('thurrun2.dat')]; % thur run 2

matdat = {dat1,dat2,dat3,dat4,dat5,dat6,dat7,dat8}; % cell for loop
% preallocating loop variables
pdat = (1:length(dat1))';
idat = pdat;
% loop to combine data
for i = 1:8 % 8 data files
    file = matdat{i}; % current element in matdat
    pdat(:,i) = file(:,1); % column of pressure
    idat(:,i) = file(:,2); % column of current
end
% averaging data
bigdat = [mean(pdat,2),mean(idat,2)]; % averages the rows of the pressure and current and then concatenates them


%p0 = 101325; % SSL pressure, Pa
p0 = 97773;
rho = 1.185; % Given density, kg/m^2
visc = 1.831e-5; % air dynamic viscosity
d = 0.2032; % given diameter of sphere (acts as chord), meters

p = bigdat(:,1).*47.88; % WT pressure
I = bigdat(:,2); % amperage from transducer used to get pressure in WT , mA

% used increasing velocity linear fit from last lab (ampfiti) I = 0.0123*p + 4.0695
% to get pressure total
pfit = @(x) 1./0.0123.*(x - 4.0695); % anonymous function used for repeat
deltaP = pfit(I); % pressure difference over sphere

% Calculating q and Re
% p1 is freestream dynamic pressure (q1)
v = sqrt(2.*p./rho); % velocity for Re calc
Re = rho.*v.*d./visc; % Reynolds number

% Re calcs (out of order a bit, used for plot below)
wtRe = interp1([Cp(13),Cp(14)],[Re(13),Re(14)],1.22,'linear'); % Wind tunnel Reynold's number from  
                                                               % interpolation when Cp = 1.22
critRe = 3.85e5; % Predefined experimental critical Reynold's number for sphere
TF = critRe./wtRe; % Turbulence factor for the wind tunnel
% Plotting deltaP/q vs. Re
figure(1)
Cp = deltaP./p;
plot(Re,Cp);
title('Average C_P vs. Re');
xlabel('Reynolds Number');
ylabel('Pressure Coefficient');
hold on
plot([0,wtRe],[1.22 1.22],'--r')
plot([wtRe, wtRe],[0, 1.22],'--k'); % plots vertical line at intersection point
legend('Average','C_P = 1.22','Re value = 3.0374e5')



% From Fig.3: Data for TF Percentage used are the two points closest to TF
% value of 1.2675, which is the result of the above TF calculation
x1 = 1; % starting point TF
x2 = 1.5; % ending point TF
per1 = 0; % starting point turbulence percent
per2 = 0.6; % ending point turbulence percent
percentTF = interp1([x1,x2],[per1,per2],TF,'linear'); % interpolating between the points at x-value 'TF'



% Individual plots
figure(2)
%1
p1 = dat1(:,1).*47.88;
I1 = dat1(:,2);
deltaP1 = pfit(I1);
plot(Re,deltaP1./p1);
hold on
%2
p2 = dat2(:,1).*47.88;
I2 = dat2(:,2);
deltaP2 = pfit(I2);
plot(Re,deltaP2./p2)
%3
p3 = dat3(:,1).*47.88;
I3 = dat3(:,2);
deltaP2 = pfit(I3);
plot(Re,deltaP2./p2);
%4
p4 = dat4(:,1).*47.88;
I4 = dat4(:,2);
deltaP4 = pfit(I4);
plot(Re,deltaP4./p4);
%5
p5 = dat5(:,1).*47.88;
I5 = dat5(:,2);
deltaP5 = pfit(I5);
plot(Re,deltaP5./p5);
%6
p6 = dat6(:,1).*47.88;
I6 = dat6(:,2);
deltaP6 = pfit(I6);
plot(Re,deltaP6./p6);
%7
p7 = dat7(:,1).*47.88;
I7 = dat7(:,2);
deltaP7 = pfit(I7);
plot(Re,deltaP7./p7);
%8
p8 = dat8(:,1).*47.88;
I8 = dat8(:,2);
deltaP8 = pfit(I8);
plot(Re,deltaP8./p8);

title('Individual C_P vs. Re');
xlabel('Reynolds Number');
ylabel('Pressure Coefficient');
legend('Tues Run 1','Tues Run 2','Wed1 Run 1','Wed1 Run 2','Wed2 Run 1','Wed2 Run 2','Thurs Run 1','Thurs Run 2');


