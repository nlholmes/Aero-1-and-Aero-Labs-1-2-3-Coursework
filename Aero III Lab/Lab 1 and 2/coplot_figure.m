clear all; close all; clc

% Create BS points within the figure's limits
% (don't try to plot something at x=2^30 when the xlim is at x=5...)
[X, Y] = pol2cart([0:0.01:(2.*pi)], cos([0:0.01:(2.*pi)])+cos(2.*[0:0.01:(2.*pi)]));
X = X+abs(min(X));
Y = Y+abs(min(Y));
X = X.*(10^6-10^5)./(max(X)-min(X));
Y = Y.*0.5./(max(Y)-min(Y));
% Less idiotic data for you fools out there:
%X = [0, linspace(10^5, 8*10^5, 20), 10^6];
%Y = [0, 0.1+0.3.*rand(1, 20), 0.5];

%% MLG 3000 patent pending code:
% Use the code from here, with your own X and Y data, saved as 'X' and 'Y'
p_X = [107, 890]; % Image's X limits, in pixels, of the plot's position
p_Y = [88, 746]; % Same for Y
X_lim = [10^4, 10^6]; % X axis limits on the figure you want to plot
Y_lim = [0, 0.5];

X_plot = p_X(1) + X.*diff(p_X)./diff(X_lim); % Magic to scale things onto plot
Y_plot = p_Y(1) + Y.*diff(p_Y)./diff(Y_lim);

figure('Name', 'MLG 3000 patent pending code', 'NumberTitle', 'off');
Fig = flip(imread('MLG3000_plot.JPG'), 1);
[stdby1, stdby2, ~] = size(Fig);
imshow(Fig);
set(gca, 'YDir', 'normal')
line(X_plot, Y_plot, 'color', 'r', 'lineWidth', 2);
print(gcf, 'MLG3000_output.png', '-dpng', '-r300');