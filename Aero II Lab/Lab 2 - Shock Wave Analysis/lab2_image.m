%% Aero II Lab, Lab 2: Image processing for Mach number
%% Load all files
% Get file name as string
% Get block number and store as array of values (for mach number calculation)
% Also get theta number and do the same
% Read each image
numFolders = 35;
numPhotos = 20;

% From Lab 2 Psuedocode and then edited
% Find all files in wdir matching the pattern DSC_*.jpg
blockNumber = [];
theta = [];
currDir = dir; % gets current directory info
for k = 1:numFolders % number of folders in the data set, skips the '.' and '..' folders
    % Get file name
    folderName{k} = currDir(k+2).name; % gets folder names and makes cell array
                                       % k+2 for the '.' and '..' folders again
    % Get block number and theta values
    blockNumber{k} = folderName{k}(6:9); % location of block number in folder name strings
    theta{k} = folderName{k}(16:end);% location of theta angle in folder name strings
    images = dir(sprintf('%s/DSC_*.jpg', folderName{k}));
    
    
    % Load all images
    for j = 1 : numPhotos % number of pictures per set
        % file name changes for each picture
        fileName = [folderName{k}, '/' , images(j).name]; % string concatenation for file path name for imread
        pics{k}{j} = imread(fileName); % loads the current image and stores in into pics cell array
                                       % k separates block number
                                       % j index for picture number
    end
end
% Convert block number and theta character cell vector values to normal arrays of numbers
blockNumber = str2num(cell2mat(blockNumber'));
theta = str2num(cell2mat(theta'));

%% Edge detection and shock line calculations
trim_val = .13; % make like .14 for this image
for k = 1:numFolders
    for j = 1:numPhotos % number of images per set, excluding first and last 6
        % Edge Detection
        % From Lab 2 Psuedocode
        pic = pics{k}{j}; % variable for reference to current image
        variation1 = diff(pic);
        variation2 = diff(pic, 2);
        edgeTol1 = max(mean(variation1));
        edgeTol2 = max(mean(variation2));
        threshold = mean([edgeTol1, edgeTol2])*trim_val; % change the .18 to some value (trim_val)
        
        % Perhaps an alternate method to line gaps -- skip bad photos
        %{
        if sum(variation1 >= 10) ~= 0 %|| sum(variation2 >= 10) ~= 0
            continue % skip this photo if there is a gap in the shock line
        end
        %}
        
        pic_cannys{k}{j} = edge(pic, 'canny', threshold); % detect edges using canny method
        
        % Cropping
        pic_cannys{k}{j} = pic_cannys{k}{j}(100:end,:); % ignore top 250 rows of pixels
        pic_canny = pic_cannys{k}{j}; % varaible for reference to current edge detected image
        %figure(2)
        %imshow(pic_canny);
        
        % Get equation for shock lines
        % Find the uppermost white pixel in each column
        %   pic_canny is a logical matrix representing the pic, with 1 where there is white
        [rows,cols] = size(pic_canny);
        colDex = []; % col of first white pixel index; array for holding the indicies of the first white pixels
        rowDex = []; % row of first white pixel index; want both rows and cols for linear fit
        count = 1; % preallocate count variable for firstDex array
        dist = 0; % preallocate distance threshold value
        colDex(1) = 0; % preallocate
        rowDex(1) = 0; % preallocate
        
        %limited_range = cols/5:1:(cols - cols/5);
        for i = floor(cols/5):1:(cols - floor(cols/5)) % for every column
            if sum(pic_canny(:,i) ~= 0) % only update array if there's a white pixel in this col
                % get first white point values
                colTemp = i;
                rowTemp = rows - find(pic_canny(:,i), 1, 'first'); % find the first nonzero (starts at top of image)
                                                                   % subtract this value from rows for plotting reasons

                % updating dist
                % only use a point that is within 10 pixels distance of the previous point used
                if count ~= 1 % do this for all except first one found
                    %colDex(count - 1)
                    %rowDex(count - 1)
                    dist = sqrt((colTemp - colDex(count - 1))^2 + ...
                        (rowTemp - rowDex(count - 1))^2);
                end

                if dist <= 10 % only update arrays if the point proximity condition is satisfied
                    colDex(count) = colTemp;
                    rowDex(count) = rowTemp;

                    count = count + 1; % update count for next index
                end
            end
            % Shock fitting and values
            shockFit{k}{j} = polyfit(colDex,rowDex,1);
            shockVals{k}{j} = polyval(shockFit{k}{j},colDex);
        end
        
        
    end
end

% White clusters will appear away from shock depending on trim_value (lower
%   value means more clusters will get through, cropping will help this
%   issue)
% Cropping to predefined size
%   Find where the white shock lines are and adda buffer to top and bottom,
%   then remove all other image material beyond the top and bottom buffer



%imshow(pic_canny);
%% Calculate shock angle from equation
%{
%   Use equation from last lab to get the Mach number from the block number
%   For the image cycling loop extract the block number which is used for
%   the shock angle calculation
%theta = acot(tan(beta)*((gamma+1)*M^2/(2*(M^2*sin(beta) - 1)) - 1);
% theta is wedge angle, we want beta for shock angle

% attempt with the equation
gamma = 1.4; % specific heat ratio of air for T < 1000K
%shockAngle = 1:length(blockNumber);
for i = 1:length(blockNumber) % calculate mach number for every block number
    M(i) = 1.82e-7 * blockNumber(i)^2 - 1.3e-3 * blockNumber(i) + 3.9; % Isentropic approximation equation for NCSU supson WT
                                                                         % BN is block number
    syms beta
    eqn = cotd(theta(i)) == tand(beta)*((gamma+1)*M(i)^2/(2*(M(i)^2*sind(beta) - 1)) - 1);
    currBeta = double(vpasolve(eqn,beta));
    if isempty(currBeta) == 1
        currBeta = 0;
    end
    shockAngle(i) = currBeta;
end
%}
%% attempt with the fit line -- this seems to give more reasonable values
%betaLine = atand(shockVals{end}{end}(end)./colDex(end));
betaLine = {};
bestPics = {};
betaVariance = {};
betaDiff = [];
%bestBeta = [];
for k = 1:numFolders % number folders (mach and theta number)
    for j = 1:numPhotos % number images per folder, but exclude the first and last 6
        % All points will have same slope, only need to get angle once
        %xdist = 1:length(shockVals{k}{j}(end));
        %betaLine{k} = atand(shockVals{k}{j}(end)/xdist(end)); % will only have a beta val for each mach number
        
        for i = 1:length(shockVals{k}{j}) % for every point on line
            %xdist = 1:length(shockVals{k}{j}(i)); % number of pixel columns
            %betaLine{k}{j}(i) = atand(shockVals{k}{j}(i)/xdist);
            betaLine{k}{j} = abs(atand(shockFit{k}{j}(1))); % gets slope from each fit
            meanBeta = mean(cell2mat(betaLine{k})); % mean of current beta set values
            %betaVariance{k} = diff(cell2mat(betaLine{k}));
            betaDiff = abs(cell2mat(betaLine{k}) - meanBeta); % difference between mean and current value of beta
            minBetaDiff = min(betaDiff); % minimum difference in beta
            %bestBeta{k} = betaLine{k}{betaDiff == minBetaDiff(k)};
            
        end
    end
    % Calculating best picture to use
    bestBeta(k) = betaLine{k}{10};
    %bestPics{k} = min(betaLine{k}{:});
end

%% Mach number calculation
tbMach = [];
for i = 1:length(blockNumber)
    M(i) = 1.82e-7 * blockNumber(i)^2 - 1.3e-3 * blockNumber(i) + 3.9; % Isentropic approximation equation for NCSU supson WT
                                                                      % BN is block number
    
    % theta beta mach relation method
    gamma = 1.4;
    syms tbM
    eqn = cotd(theta(i)) == tand(bestBeta(i))*((gamma+1).*tbM.^2./(2.*(tbM.^2.*sind(bestBeta(i)) - 1)) - 1);
    mSol = double(vpasolve(eqn,tbM));
    mSol = mSol(mSol > 0); % take only positive solution
    tbMach{i} = mSol;
end
M = double(M);

tbMach = tbMach(~cellfun('isempty',tbMach)); % remove empty cells
tbMachAvg = [];
p = 1;
for i = 1:7 % number of wedge angles times block number (mach number)
    tbMachAvg(i) = mean(cell2mat(tbMach(p:p+3))); % only have 4 mach numbers per angle now
    p = p + 4; % update index variable
end
%%
figure(1)
plot(theta(1:5),bestBeta(1:5),'*-')
hold on
plot(theta(6:10),bestBeta(6:10),'*-')
plot(theta(11:15),bestBeta(11:15),'*-')
plot(theta(16:20),bestBeta(16:20),'*-')
plot(theta(21:25),bestBeta(21:25),'*-')
plot(theta(26:30),bestBeta(26:30),'*-')
plot(theta(31:35),bestBeta(31:35),'*-')
mach = {};
for i = 1:7
    mach{i} = sprintf('M = %.3s',M(i*5));
end
title('\theta-\beta-M Oblique Shock Angle Relation')
legend(mach{1:7}); % data is block number separated (mach number representation)
xlabel('Wedge Angle (\theta, deg)');
ylabel('Shock Angle (\beta, deg)');
%{
%for f = 1:7 % number of unique block numbers
   while q <= 31 % number of unique wedge angles
       plot(theta(q:q+4), bestBeta(q,q+4), '*-');
       hold on
       q = q + 1;
   end
   hold off
%end
%}
%plot(theta,-shockAngle)
hold off
%% plot mach number vs block number
tbBlockNumber = [];
for i = 1:5:35
    tbBlockNumber = [tbBlockNumber; blockNumber(i+1:i+4)];
end
figure(2)
plot(tbBlockNumber,cell2mat(tbMach)) % discontinuities due to solver clearly seen, mach lower
hold on
plot(blockNumber,M);
xlabel('Block Number');
ylabel('Mach Number'); 
title('Mach Number vs. Block Number');
legend('\theta-\beta-M relation','Lab 1 calibration equation');
