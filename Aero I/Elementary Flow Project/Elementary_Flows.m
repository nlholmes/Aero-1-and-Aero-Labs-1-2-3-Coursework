%% Needed Corrections
% Change all the lims and stuff to be multed by v_inf so higher vels work
% Get the section lines in the flood plot to disappear so just color
% Add colorbar to flood plot as a legend
% Fix lifting cylinder stagnation streamline
% Fix velocity plots not working properly (could be the v_inf multiply issue as that fixed psi)
%% Coordinates
% Set Bounds limit
x_lim = 2;
y_lim = 0.5*x_lim;

v_inf = 8; % V_inf moved up here as broken-ness

% Non-cylinder XY space (starts from 0)
x0 = 0;
y0 = 0;
x = linspace(0,x_lim,200); 
y = linspace(0,y_lim,100); 
[X,Y] = meshgrid(x-x0,y-y0);

% Cylinder XY space (starts from -xlim)
xcyl = linspace(-x_lim,x_lim,200); % x lim: make -2 to 2 for cylinder 
ycyl = linspace(-y_lim,y_lim,100); % y lim: make -1 to 1 for cylinder 
[Xcyl, Ycyl] = meshgrid(xcyl,ycyl);

%% Input Variables
              % Freestream velocity
dist = 1;               % distance between source/sink
lambda = 1*v_inf;             % Volumetric Flow
k = lambda*dist;        % strength factor

% Cylinder specific
epsilon = 1e-6; % Adjustment factor to prevent r, rcyl from going to zero
%R = 0.5;
R = sqrt(k./(2*pi*v_inf));  % used for simplification
rcyl = sqrt((Xcyl.^2 + Ycyl.^2 + epsilon.^2)); % current radius of cylinder
% Lifting Cylinder Only
constant = 0.5;             % Constant for gamma value, increasing this will lower the positions of the stagnation points
                            % Zero results in a regular non-lifting cylinder
gamma = 2*pi*constant;      % Circulation value, for use in point vortex

% Source Position
xsource = 0.5;          % source position is inputted
ysource = 0.5;          
% Sink Position
xsink = xsource + dist; % sink position is the specified distance away
ysink = 0.5;            % ysource = ysink at 0 alpha gives oval

% Velocity Calculations
rsource = sqrt(((X - xsource).^2 + (Y - ysource).^2 + epsilon.^2)); % r^2 = x^2 + y^2 for velocity conversion to cartesian and calculations
rsink = sqrt(((X - xsink).^2 + (Y - ysink).^2 + epsilon.^2));

%% Uniform Flow
psi_unif = v_inf*(Y-ysource);               % Y - ysource gives the correct psi = 0
u_unif = v_inf;                             % From notes
v_unif = 0;
vel_unif = sqrt(u_unif.^2 + v_unif.^2);     % For Velocity Plots
%% Source
psi_source = lambda./(2*pi).*atan2(Y - ysource, X - xsource);   % adjusts for source position
u_source = lambda./(2*pi).*X./(rsource.^2);
v_source = lambda./(2*pi).*Y./(rsource.^2);
vel_source = sqrt(u_source.^2 + v_source.^2);
%% Sink
psi_sink = -lambda./(2*pi).*atan2(Y - ysink, X - xsink);        % adjusts for sink position
u_sink = -lambda./(2*pi).*X./(rsink.^2);
v_sink = -lambda./(2*pi).*Y./(rsink.^2);
vel_sink = sqrt(u_sink.^2 + v_sink.^2);
%% Rankine Half-Oval (unif + source)
psi_halfoval = psi_unif + psi_source;
u_halfoval = u_unif + u_source;
v_halfoval = v_unif + v_source;
vel_halfoval = sqrt(u_halfoval.^2 + v_halfoval.^2);
%% Rankine Oval (unif + source + sink)
psi_oval = psi_unif + psi_source + psi_sink;    % need subtract ysource for the vertical position?
                                                % it works but is it what I am supposed to do?
                                                % Solved?: instead subtracted ysource from unif flow Y value
u_oval = u_unif + u_source + u_sink;
v_oval = v_unif + v_source + v_sink;
vel_oval = sqrt(u_oval.^2 + v_oval.^2);
%% Doublet
psi_doublet = -k./(2*pi).*Ycyl./rcyl.^2;
u_doublet = v_inf.*dist.*(Ycyl.^2 - Xcyl.^2 - epsilon^2)./rcyl.^4; % dist = l
v_doublet = -v_inf.*dist.*2.*Xcyl.*Ycyl./rcyl.^4;
vel_doublet = sqrt(u_doublet.^2 + v_doublet.^2);
%% Point Vortex
psi_vort = gamma./(2*pi).*log(rcyl); % only used for cylinders
u_vort = gamma./(2*pi).*Ycyl./(rcyl.^2);
v_vort = -gamma./(2*pi).*Xcyl./(rcyl.^2);
vel_vort = sqrt(u_vort.^2 + v_vort.^2);
%% Non-Lifting Cylinder (unif + doublet)
% Cylinders require xlim and ylim changes (top of script)
psi_cyl = psi_unif + psi_doublet;

u_cyl = v_inf.*(1 + R.^2*(rcyl.^2 - 2.*Ycyl.^2)./rcyl.^4); % Worked out on paper first, then adjusted for variables
v_cyl = -v_inf.*R.^2.*2.*Xcyl.*Ycyl./rcyl.^4;
vel_cyl = sqrt(u_cyl.^2 + v_cyl.^2);
%% Lifting Cylinder (unif + doublet + vort)
%psi_liftcyl = psi_cyl + gamma/(2*pi)*log(rcyl./R);
psi_liftcyl = psi_cyl + psi_vort;

u_liftcyl = u_cyl + u_vort;
v_liftcyl = v_cyl + v_vort;
vel_liftcyl = sqrt(u_liftcyl.^2 + v_liftcyl.^2);
%% 'Invented' Flows
psi_invent1 = psi_unif + psi_source + psi_sink + psi_doublet + psi_vort; % v_inf = 0.25, dist = 0.1, gamma const. = 0.9
u_1 = u_liftcyl + u_source + u_sink;
v_1 = v_liftcyl + v_source + v_sink;
vel_invent1 = sqrt(u_1.^2 + v_1.^2);
psi_invent2 = psi_unif + psi_source + psi_sink + psi_doublet; % lambda = 0.2, dist = 0.65,
                                                              % xsource = 0.9, ysource = 0.3, ysink = 0.1
                                                              % gamma const = 0.4, v_inf = 0.25
u_2 = u_cyl + u_source + u_sink;
v_2 = v_cyl + v_source + v_sink;
vel_invent2 = sqrt(u_2.^2 + v_2.^2);
%% Plotting (choose which to plot)
% Stream Function (psi)
% Which to plot?
selectedPsi = psi_oval;
selectedVel = vel_oval;
stagpt = 0; % zero unless half oval, which it is then lambda/2, represents stagnation streamline
% Stream Function
psiStep = 0.1*v_inf; % step for which contours to plot, 0.1 good for non-cyl, 0.3 good for cyl
               % can remove labels to reduce cluster
velStep = 0.1; % step for velocity plot
figure(1)
psiplot = contour(X,Y,selectedPsi,[-3*v_inf:psiStep:3*v_inf],'ShowText','on','LineColor','k');
pbaspect([2 1 1]); % correct aspect ratio [x y z]
hold on
stagplot = contour(X,Y,selectedPsi,[stagpt,stagpt],'LineColor','r','LineWidth',1.5); % shows the 'body' where psi = 0
stagplot2 = contour(X,Y,selectedPsi,[-stagpt,-stagpt],'LineColor','r','LineWidth',1.5); % useful only when half oval
hold off

% Velocity Plot (flood plot) (vel)
figure(2)
velplot = contourf(X,Y,selectedVel,[-3:velStep:3],'showtext','on','labelspacing',288);
hold on
stagplotvel = contour(X,Y,selectedPsi,[stagpt,stagpt],'LineColor','r','LineWidth',1.5); % based off psi, shows body
stagplotvel2 = contour(X,Y,selectedPsi,[-stagpt,-stagpt],'LineColor','r','LineWidth',1.5); % based off psi, shows body
pbaspect([2 1 1]);
hold off


