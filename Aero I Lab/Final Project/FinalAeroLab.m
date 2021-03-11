%% Final Aero I Lab Project: Drag stuff my guy
% y pos is up
% z pos is how far back wake rake is
% each z_loc file has 6 y changes up to 2 ft total change, with static zloc
% files are in inches, psf, english units

% Loop will be required for each zloc file for the 6 yloc readings smooshed
% Loop will also be required to get the file name --> do nested loops? nah
%   Only concern with file name for first: the (1) etc
%   Second loop for the yloc changes: do separate loop to have big files
%       for error, then separate into average files with each of the 10
%       readings per yloc averaged
% Worry about averages and stuff after the files are separated to speed up
%   the code (use matrices)


% Could just use load as welladadaddadad
%% Code
% Adding file paths
folder = fileparts(which('FinalAeroLab.m'));
addpath(genpath(folder));

% Loop for big file variables (textscan)

numDynP = 6; % there are 6 dynamic pressures that data was gathered for
zLocs = 0:15; % z locations as numbered in each file(0 thru 15)
formatspecs = '';
lines = 70; % number of columns in each
for w = 1:lines % for number of lines in file
    formatspecs = [formatspecs,'%f']; % makes the delimiter
end
bigfileset = {}; % preallocates bigfile
for i = 1:numDynP % repeat for each dynamic pressure data was gathered at
    for j = zLocs
        % Condition for file naming, the first does not have a (#)
        if i == 1 % first folder different
            file = ['z-loc_',sprintf('%d',j),'.txt']; % gives current file name in .txt
        else 
            file = ['z-loc_',sprintf('%d',j),sprintf(' (%d)',i - 1),'.txt']; % adds the (#), starts at j - 1
        end
    end
    fid = fopen(file);
    bigfileset{i} = textscan(fid,formatspecs);
    % want file so that the rows are the q's and in them are the 70 cells
    bigfile = 0; % corrected version
    
end
% Loop for separate ylocs (take certain elements from big file variables)
