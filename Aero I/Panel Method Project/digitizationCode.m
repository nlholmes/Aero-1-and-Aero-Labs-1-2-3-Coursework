clear all; close all; clc
% Anton Scholten, NCSU
% This basic code will digitalize the contour of an image and give you the
% x, y coordinates of this contour.
% A file named 'Output_edge_coordinates.txt' will be created where you run
% this code and have the number of panels, x-y coordinates.
% IF YOUR PICTURE IS WEIRD, the results may be weird

% For best results, import a JPG picture that has:
% _ white (or light in color) background
% _ no background color within the body of the shape to detect
% _ dark edges
% _ square aspect ratio
% Ideal image is something similar to a clipart

pic = imread('Squid_jet_outline_gimp_fill.jpg'); % Change the name to your picture
threshold = 245; % Value between 0 (black) and 255 (white) below which grayscale values are set to zero (black)
num_coordinates = 200; % How many points to save

figure('Name', '', 'NumberTitle', 'off', 'units', 'normalized', 'outerposition', [0.2, 0.2, 0.6, 0.7])
subplot(2, 3, 1)
imshow(pic) % Picture as it is
title('Original')

pic2 = rgb2gray(pic);
subplot(2, 3, 2)
imshow(pic2) % Grayscale of the picture
title('Grayscale')

% Fill grayscale
subplot(2, 3, 3)
pic3 = pic2;
pic3(pic3 < threshold) = 0;
imshow(pic3) % Filled grayscale
title('Filled grayscale')

% Edge detect
pic4 = edge(pic3);
subplot(2, 3, 4)
imshow(~pic4) % Edge
title('Edge detection')

[x, y] = find(rot90(pic4, -1) == 1); % Find indices of edge
x = x./size(pic4, 2); % Normalize
y = y./size(pic4, 1);

%% Sort points based on angle
cx = mean(x); % Find centroid of shape
cy = mean(y);
a = atan2(y - cy, x - cx); % Find angle between any point and centroid
[sorted_angles, order] = sort(a); % Sort angles in ascending order and remember their index
x_sort = x(order); % Sorted values
y_sort = y(order);

% Sort by trailing edge
[~, indx] = max(x_sort);
x_sort = flipud([x_sort(indx:end); x_sort(1:indx-1)]);
y_sort = flipud([y_sort(indx:end); y_sort(1:indx-1)]);

% Bring the shape down to the x-axis
[~, indx] = min(x_sort);
y_sort = y_sort-y_sort(indx);

%% Second sorting using closest distance
% This cell of the code was inspired from:
% https://stackoverflow.com/questions/11631934/sort-coordinates-points-in-matlab
data = [x_sort, y_sort];
dist = pdist2(data, data); % Find distance between each point and the next
dist(dist > 0.01) = Inf; % Get rid of outliers

result = zeros(1, length(x_sort)); % Pre-allocate
result(1) = 1; % first point is first row in data matrix

for i=[2:length(x_sort)] % Check all the distances
    dist(:, result(i-1)) = Inf; % Set previous distances to infinity so they do not bother the next point's verification
    [~, idx] = min(dist(result(i-1), :)); % Find the next closest point and save it's index
    result(i) = idx;
end

x_sort = x_sort(result);
y_sort = y_sort(result);


%% Write onto output file
save_index = round(linspace(1, length(x_sort), num_coordinates));
fid = fopen('Output_edge_coordinates.txt', 'wt');
fprintf(fid, sprintf('%d Number of panels \n', num_coordinates-1));
fclose(fid);
dlmwrite('Output_edge_coordinates.txt', [x_sort(save_index), y_sort(save_index)], '-append', 'delimiter', '\t', 'newline', 'pc', 'precision', 4)

%% Plot to verify order of coordinates
subplot(2, 3, 5) % Plot the unsorted coordinates
p1 = plot(x(1), y(1), 'k.');
hold on
axis([-0.1, 1.1, -0.1, 1.1])
title('Raw edge coordinates')

subplot(2, 3, 6) % Plot coordinates after dual sort
p2 = plot(x_sort(1), y_sort(1), 'b.');
hold on
axis([-0.1, 1.1, -0.6, 0.6])
title('Angle and distance sorting')
tic
for i=round(linspace(1, length(x), 100)) % 300 iterations with a pause of 0.01 second = 3 seconds
    set([p1; p2], {'xdata'}, {x(1:i); x_sort(1:i)}, {'ydata'}, {y(1:i); y_sort(1:i)})
    pause(0.01)
end
hold off

figure('Name', 'Saved Coordinates', 'NumberTitle', 'off')
axis([-0.1, 1.1, -0.6, 0.6])
hold on
for i=[1:length(save_index)]
    plot(x_sort(save_index(i)), y_sort(save_index(i)), '+k')
    pause(0.01)
end
hold off
