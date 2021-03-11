%% Final Aero I Project: Panel Method over 2 airfoils and an animal
% Digitization software used for animal

% Needs:
% Coordinates for airfoil
%   files on website start backwards (do in opposite order or remove first name bit and flip the mat)
% Equations from class
% Try to do without loops --> use vectorization and matrices (later)

% New needs
% Plot the contour bit with the flow velocities (extra credit)
% Camber of living thing
% Digitize living thing
% Run throughs for 2 airfoils and a living thing
% PLOT AT THE MIDPOINTS OF THE PANELS*** and change loop variables to go off of midpoints
% cl error percentages
clear all;close all;
%% Start of code
% Loading file and adjusting
%   Remove name
%   flipud
%   Only works for semi-sign-consistent files (see flipud section, first 8 elements)

% ---------------- Change file here ---------------- %
fid = fopen('Output_edge_coordinates.txt'); % filename
% ---------------- ---------------- ---------------- %

fgetl(fid); % skip first line
rawFoil = textscan(fid,'%f %f'); % textscan it in
fclose('all');

% Flipud to get TE to LE (ONLY DO THIS IF IT IS FROM THE WEBSITE)
if sum(rawFoil{2}(1:8) < 0) <= 0          % if there is not a negative number in the first 8 elements
                                         %      (means starts bottom side)
                                         %  then it needs to be flipped upside down for panel
                                         %  method to work
    x = flipud(rawFoil{1}); % x coords are first col
    y = flipud(rawFoil{2}); % y coords are second
else % if not (animal, me163)
    x = rawFoil{1};
    y = rawFoil{2};
end

% inputs:      uinf (free-stream velocity)
%              alpha (angle of attack; degrees)
      aalphad = -5:15; % angle of attack in deg(-5 to 15 for the assignment)
      aalpha = aalphad.*pi./180.0; % degrees to radians
      uinf = 1;
      c = 1; % chord for coeff lift calculation
      npanel = length(x) - 1; % number of panels (210 for me163, each coord vec is 211 long)

%% jacks code
%pi = 4.0*atan(1.0);


figure(1) % cp plot
%c ----- import panel coordinates (x, y (1:npanel+1))
%c ----- Note:  ordering must be from bottom trailing edge to top trailing edge (clockwise)
for q = 1:length(aalpha) % for the angles of attack
    aalphaLoop = aalpha(q); % current angle of attack indexed
      for j=1:npanel
       ds(j) = sqrt((x(j+1)-x(j))^2 + (y(j+1)-y(j))^2); % panel length
       tnx(j) = (x(j+1)-x(j))/ds(j);  % x component of panel tangent = cos(theta_j)
       tny(j) = (y(j+1)-y(j))/ds(j);  % y component of panel tangent = sin(theta_j)
       xnx(j) = -tny(j);              % x component of panel normal
       xny(j) =  tnx(j);              % y component of panel normal
      end

%c ---- apply V dot n = 0.0 for every panel

      for i=1:npanel
        xi = 0.5*(x(i)+x(i+1));
        yi = 0.5*(y(i)+y(i+1));
        sumn = 0.0;
        sumt = 0.0;
        for j=1:npanel
         xj = x(j);
         yj = y(j);
         xip =  tnx(j)*(xi-xj) + tny(j)*(yi-yj);  %x* location in panel coord. system          
         yip = -tny(j)*(xi-xj) + tnx(j)*(yi-yj);  %y* location in panel coord. system
         upv = 0.5/pi*(atan2(yip,xip-ds(j))-atan2(yip,xip)); %x* velocity in panel coord. system.
         vpv = 0.25/pi*log(((xip-ds(j))^2 + yip^2)/(xip^2 + yip^2)); %y* velocity in panel coord. system
         if (i==j) 
         upv = 0.5;
         vpv = 0.0;
         end
         uv = tnx(j)*upv - tny(j)*vpv;  %x component of induced velocity in Cart. system
         vv = tny(j)*upv + tnx(j)*vpv;  %y component of induced velocity in Cart. system
         us = -vv; % x component of source velocity
         vs =  uv; % y component of source velocity
         a(i,j)  = us*xnx(i) + vs*xny(i); %matrix elements
         at(i,j) = us*tnx(i) + vs*tny(i); %matrix elements storing tangential components
         sumn = sumn + uv*xnx(i) + vv*xny(i);
         sumt = sumt + uv*tnx(i) + vv*tny(i);
        end
        
        a(i,npanel+1) = sumn;
        b(i) = -uinf*(cos(aalphaLoop)*xnx(i) + sin(aalphaLoop)*xny(i));
        at(i,npanel+1) = sumt; 
      end

%c --- apply Kutta condition 
      for j=1:npanel+1
       a(npanel+1,j) = at(1,j)+at(npanel,j);
      end % idk if b is supposed to be in this for loop or not, as linsolve reqs same size and it fixes it
      b(npanel+1) = -uinf*(cos(aalphaLoop)*(tnx(1)+tnx(npanel)) + sin(aalphaLoop)*(tny(1)+tny(npanel)));

%c ----now solve A*ss = b to get the source strengths (ss(1:npanel)) and vortex strength (ss(npanel+1))
%c     Note that your matrix is (npanel+1,npanel+1)
ss = linsolve(a,b'); % wrote this bit, had to transpose b for it to properly work

%c ---- now compute tangential velocity and cp for each panel
      for i=1:npanel
       xi = 0.5*(x(i)+x(i+1));
       yi = 0.5*(y(i)+y(i+1));
       jacksum = 0.0;
       for j=1:npanel+1
        jacksum = jacksum + at(i,j)*ss(j);
       end
       vtan = jacksum + uinf*(cos(aalphaLoop)*tnx(i) + sin(aalphaLoop)*tny(i));
       cp(i) = 1.0 - vtan^2/uinf^2; %Cp
      end
      
      
% Coefficient of Lift Calculation
% Panel Method with ss
    cl(q) = 2.*ss(end)./(c.*uinf).*sum(ds); %will have 1 cl for each alpha

% Cp plots
plot(x(1:npanel),cp');  % have to change teh axis labels as they are
                        % actually negative even though says positive
                        % exclude last coord in x
set(gca,'YDir','reverse'); % changed axis instead of changing values
hold on
end

%% TAFT cl calculation
% Need camber first
% camber = 1/2[(y/c)|up + (y/c)|low]
%   add because half of the y values are negative
%   need y/c divided into upper and lower coords, x/c stays same

% For uneven panel numbers
if rem(length(y),2) ~= 0 % if odd
    ylow = y(1:floor(1/2*length(y)));
    yup = y(ceil(1/2*length(y):end-1)); % leave out last coord
else % if even
    ylow = y(1:1/2*length(y));
    yup = y(1/2*length(y)+1:end);
end

ycamb = 1/2.*(yup + ylow);
% theta calculation
% theta = acos[1-2x(i)] % transformation eqns
% theta does not suffer the plotting problem, see later comments on the airfoil plot
xcamb = x(ismember(y,yup));
theta = acos(1 - 2.*xcamb);

% Actual cl calc
% cl = 2*pi[alpha - 1/pi*sum(df/dx(cos(theta) - 1)*deltaTheta]
% df/dx ~~ deltaY/deltaX
% cos(theta) - 1 = -2*1/2(x(i) - x(i+1)), where the 1/2 onward part is the x midpt of the panel
% Jack's method: numerical loop
%   could try trapz or syms int instead
for q = 1:length(aalpha)
    alphaLoop = aalpha(q); % same setup as before for current alpha
    for i = 1:(length(yup)-1) % number of panels (number of midpoints)
        dfdx = (ycamb(i+1) - ycamb(i))./(xcamb(i+1) - xcamb(i));
        deltaTheta = theta(i+1) - theta(i);
        thetaBits = -(xcamb(i) + xcamb(i+1));
        bigSum = sum(dfdx.*thetaBits.*deltaTheta);
    end
    clTAFT(q) = 2*pi*(alphaLoop - 1/pi.*bigSum);
end

%% Plotting
title('C_p vs. x/c')
xlabel('x/c')
ylabel('C_p')
hold off

% Airfoil Plot
figure(2) % airfoil alone plot, with camber
plot(x,y); % airfoil
title('Airfoil')
ylim([-1,1])
xlabel('x/c')
ylabel('y/c')
hold on
plot(linspace(0,max(x),length(ycamb)),ycamb,'r--'); % camber plot
%plot(xcamb,ycamb,'r--');
hold off
%eeee = x(ismember(y,yup));
%eeee(1) = 0;
% ismember gives problem plotting because values for x start and end at 1
%plot(eeee, camb, 'k--')

figure(3) % cl plot
% Panel method
plot(aalphad,cl)
hold on
% TAFT
plot(aalphad, clTAFT)
title('C_L vs. alpha');
xlabel('Alpha (deg)')
ylabel('C_L')
legend('Panel Method','TAFT');

%% Data manipulaion
clDiff = abs(clTAFT - cl); % difference in cls
clVar = mean(clDiff); % average variation in cls

% max cl
clMax = max(cl);
clTAFTMax = max(clTAFT);
% max cp
cpMax = min(cp); % final angle of attack
%% Extra Credit
%{
% Flow at a point away from the body
% Contour plots return
% Meshgrid thing
% Equation from class
% Matlab streamline function and one other function needed
% Velocity = vinf + vinduced
x_lim = 2;
y_lim = 0.5*x_lim;
x0 = 0;
y0 = 0;
x = linspace(0,x_lim,200); 
y = linspace(0,y_lim,100); 
[X,Y] = meshgrid(x-x0,y-y0);



for q = 1:length(aalpha) % for the angles of attack
    aalphaLoop = aalpha(1); % current angle of attack indexed
    uinfx(q) = uinf*cos(aalphaLoop); % x direction freestream vel
    uinfy(q) = uinf*sin(aalphaLoop); % y direction freestram vel
      for j=1:npanel
       ds(j) = sqrt((x(j+1)-x(j))^2 + (y(j+1)-y(j))^2); % panel length
       tnx(j) = (x(j+1)-x(j))/ds(j);  % x component of panel tangent = cos(theta_j)
       tny(j) = (y(j+1)-y(j))/ds(j);  % y component of panel tangent = sin(theta_j)
       xnx(j) = -tny(j);              % x component of panel normal
       xny(j) =  tnx(j);              % y component of panel normal
      end

%c ---- apply V dot n = 0.0 for every panel

      for i=1:npanel
        xi = 0.5*(x(i)+x(i+1));
        yi = 0.5*(y(i)+y(i+1));
        sumn = 0.0;
        sumt = 0.0;
        for j=1:npanel
         xj = x(j);
         yj = y(j);
         xip =  tnx(j)*(xi-xj) + tny(j)*(yi-yj);  %x* location in panel coord. system          
         yip = -tny(j)*(xi-xj) + tnx(j)*(yi-yj);  %y* location in panel coord. system
         upv = 0.5/pi*(atan2(yip,xip-ds(j))-atan2(yip,xip)); %x* velocity in panel coord. system.
         vpv = 0.25/pi*log(((xip-ds(j))^2 + yip^2)/(xip^2 + yip^2)); %y* velocity in panel coord. system
         if (i==j) 
         upv = 0.5;
         vpv = 0.0;
         end
         uv = tnx(j)*upv - tny(j)*vpv;  %x component of induced velocity in Cart. system
         vv = tny(j)*upv + tnx(j)*vpv;  %y component of induced velocity in Cart. system
         us = -vv; % x component of source velocity
         vs =  uv; % y component of source velocity
         a(i,j)  = us*xnx(i) + vs*xny(i); %matrix elements
         at(i,j) = us*tnx(i) + vs*tny(i); %matrix elements storing tangential components
         sumn = sumn + uv*xnx(i) + vv*xny(i);
         sumt = sumt + uv*tnx(i) + vv*tny(i);
         
         vel(j,:) = sqrt((uv + uinfx).^2+(vv + uinfy).^2);
        end
        
        a(i,npanel+1) = sumn;
        b(i) = -uinf*(cos(aalphaLoop)*xnx(i) + sin(aalphaLoop)*xny(i));
        at(i,npanel+1) = sumt; 
      end

%c --- apply Kutta condition 
      for j=1:npanel+1
       a(npanel+1,j) = at(1,j)+at(npanel,j);
      end % idk if b is supposed to be in this for loop or not, as linsolve reqs same size and it fixes it
      b(npanel+1) = -uinf*(cos(aalphaLoop)*(tnx(1)+tnx(npanel)) + sin(aalphaLoop)*(tny(1)+tny(npanel)));

%c ----now solve A*ss = b to get the source strengths (ss(1:npanel)) and vortex strength (ss(npanel+1))
%c     Note that your matrix is (npanel+1,npanel+1)
ss = linsolve(a,b'); % wrote this bit, had to transpose b for it to properly work

%c ---- now compute tangential velocity and cp for each panel
      for i=1:npanel
       xi = 0.5*(x(i)+x(i+1));
       yi = 0.5*(y(i)+y(i+1));
       jacksum = 0.0;
       for j=1:npanel+1
        jacksum = jacksum + at(i,j)*ss(j);
       end
       vtan = jacksum + uinf*(cos(aalphaLoop)*tnx(i) + sin(aalphaLoop)*tny(i));
       cp(i) = 1.0 - vtan^2/uinf^2; %Cp
      end
end

contour(X,Y,vel(:,1));
%}
fclose('all');