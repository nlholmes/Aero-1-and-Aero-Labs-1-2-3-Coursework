% cylinder.m 12/1/2009 Parker MacCready
%
% this plots the streamfunction and velocity potential for potential flow
% around a cylinder
clear
% make axes
xymax = 2;
x = linspace(-xymax,xymax,100);
y = linspace(-xymax,xymax,100);
% note that x and y don't include 0
[X,Y] = meshgrid(x,y);
R = sqrt(X.^2 + Y.^2);
sin_th = Y./R;
cos_th = X./R;
U = 1;
a = 1;
phi = U*(R + a*a./R).*cos_th;
psi = U*(R - a*a./R).*sin_th;

psi1 = U.*Y.*(1-a.^2./R.^2);


figure
contour(X,Y,phi,[-3:.25:3],'-r');
hold on
[cc,hh] = contour(X,Y,phi,[-3:1:3],'-r');
clabel(cc,hh);
contour(X,Y,psi1,[-3:.25:3],'-b');
[cc,hh] = contour(X,Y,psi,[-3:1:3],'-b');
clabel(cc,hh);
xlabel('X (m)')
ylabel('Y (m)')
title('\phi=RED \psi=BLUE')
axis equal
axis tight