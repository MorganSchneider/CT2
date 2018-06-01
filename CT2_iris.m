% This script will be used to derive the temperature structure-function
% parameter from Iris+ data.
%
% Partial credit to Chilson bc I copy/pasted basically the first quarter
%
% Now that the Iris+ has been phased out in favor of the Coptersonde(s),
% this code is pretty much obsolete.  It was a good run, Iris, but on to
% bigger and better things.

%% Set up & user inputs
clear all
clc

% Enter date of flight
procYear = 2017;
procMonth = 03;
procDay = 25;

% Enter desired lat/lons
% get waypoints from mission planner - do this for coptersonde
desiredLats = [34.982715, 34.982655, 34.982490, 34.982265, 34.982035, 34.981875, 34.981815, 34.981877, 34.982040, 34.982266, 34.982490, 34.982658];
desiredLons = [-97.52109, -97.52081, -97.52061, -97.520545, -97.520615, -97.52082, -97.52109, -97.52137, -97.521565, -97.52164, -97.521565, -97.52136];
if length(desiredLats) ~= length(desiredLons)
    error('Desired lat and lon vectors are not the same length');
end
nClusters = length(desiredLats);

% Flag to decide if image file should be created
imgFlag = true;

% *** You will need to change the baseDir for your computer
% This is where your 'thermo' folder lives
baseDir = '/Users/Morgan/Documents/MATLAB/CLOUDMAP/';

%% Read in the data

% Create the directory of the matlab library and add it to the path
libDir = [ baseDir 'Github' filesep 'thermo' filesep 'matlab' filesep ];
addpath(libDir)

% Read the copter data
vehicleType = 'iris+';
% Find the appropriate directory based on instrument type
dataDirName = getDataDir(baseDir, procYear, procMonth, procDay, vehicleType);
inFileName = sprintf('%4.4d%2.2d%2.2d.mat', procYear, procMonth, procDay);
fprintf('Reading file: %s%s\n', dataDirName, inFileName)
load([ dataDirName inFileName ]);

% Read the sensor data
iFlight = input('Enter the flight number to process: ');
sensorType = 'Windsond';
% Find the appropriate directory based on instrument type
dataDirName = getDataDir(baseDir, procYear, procMonth, procDay, sensorType);
% Get the log file
logFileName = sprintf('%4.4d%2.2d%2.2d_log.txt', procYear, procMonth, procDay);
fp = fopen([ dataDirName logFileName ]);
while ~feof(fp);
    str = fgetl(fp);
    if strfind(str, '#Flight')
        logFlight = fscanf(fp, '%d', 1);
        nSensorFiles = fscanf(fp, '%d', 1);
        if logFlight == iFlight
            for iFile = 1: nSensorFiles
                sensorFileName{iFile} = fscanf(fp, '%s', 1);
            end
            break
        else
            for iFile = 1: nSensorFiles
                fgetl(fp);
            end
        end
    end
end

for iFile = 1: nSensorFiles
    [windsond, ~] = readWindsond(dataDirName, sensorFileName{iFile});
    windsondArr(iFile) = windsond;
end

nNans = 10; % arbitrarily chosen
for iFile = nSensorFiles + 1: 4
    windsondArr(iFile).obsTime = nan(1, nNans);
    windsondArr(iFile).battery_V = nan(1, nNans);
    windsondArr(iFile).pressure_Pa = nan(1, nNans);
    windsondArr(iFile).temperature_C = nan(1, nNans);
    windsondArr(iFile).humidity_perCent = nan(1, nNans);
    windsondArr(iFile).temperatureSensor_C = nan(1, nNans);
    windsondArr(iFile).latitude_deg = nan(1, nNans);
    windsondArr(iFile).longitude_deg = nan(1, nNans);
    windsondArr(iFile).altitude_m = nan(1, nNans);
    windsondArr(iFile).speed_mps = nan(1, nNans);
    windsondArr(iFile).heading_deg = nan(1, nNans);
end

lenArr = [];
for n = 1:nSensorFiles
    l = length(windsondArr(n).obsTime);
    lenArr = [lenArr l];
end
lenStr = num2str(lenArr);
fprintf('Sensor sample points: %s \n', lenStr);
prompt = 'Enter indices of sensor files to ignore, separated by spaces: ';
indStr = input(prompt,'s');
indNum = str2num(indStr);
none = isempty(indNum);
if ~none
    for m = 1:length(indNum)
        windsondArr(indNum(m)).obsTime = nan(1, nNans);
        windsondArr(indNum(m)).temperature_C = nan(1, nNans);
        windsondArr(indNum(m)).humidity_perCent = nan(1, nNans);
        windsondArr(indNum(m)).temperatureSensor_C = nan(1, nNans);
    end
end


% Remove the matlab library
rmpath(libDir)

%% Find time range. (Interactive)

timeBeg = nanmin([windsondArr(1).obsTime(1) windsondArr(2).obsTime(1) windsondArr(3).obsTime(1) windsondArr(4).obsTime(1)]);
timeEnd = nanmax([windsondArr(1).obsTime(end) windsondArr(2).obsTime(end) windsondArr(3).obsTime(end) windsondArr(4).obsTime(end)]);
indGPS = find(timeBeg <= iris.obsTimeGPS & iris.obsTimeGPS <= timeEnd);

nCircles = input('Input number of flight circles to process: ');
figure(1)
clf
plot(iris.obsTimeGPS(indGPS), iris.latitudeGPS_deg(indGPS))
xlabel('Time UTC')
ylabel('Iris GPS Latitude (deg)')
datetick('x', 13)
fprintf('Click circle start and end times\n')
shg

timeCircStart = [];
timeCircEnd = [];
for i = 1:nCircles
    [x, ~] = ginput(2);
    starts = x(1);
    ends = x(2);
    if (starts < timeBeg || starts > timeEnd), starts = nan; end
    if (ends < timeBeg || ends > timeEnd), ends = nan; end
    timeCircStart = [timeCircStart starts];
    timeCircEnd = [timeCircEnd ends];
end

%% Find mean temperature of each waypoint

% Create lat/lon clusters
latWindow = 0.00003; %these are not necessarily universal values
lonWindow = 0.00003;
latMin = desiredLats - latWindow;
lonMin = desiredLons - lonWindow;
latMax = desiredLats + latWindow;
lonMax = desiredLons + lonWindow;

% Find temps (each iteration finds cluster temps for one flight).
flight = struct('times',[],'temps',[],'latMean',[],'lonMean',[],'tempMean',[],...
    'xMean',[],'yMean',[],'radMean',[],'circumference',[],'tempDiffs',[]);

for i = 1:nCircles
    % This finds the indices of all times between the circle start and end
    % times (basically, the times during which the vehicle is flying in the
    % circle of constant height) and puts those indices into their own array.
    indWindsond1 = find(timeCircStart(i) <= windsondArr(1).obsTime & windsondArr(1).obsTime <= timeCircEnd(i));
    indWindsond2 = find(timeCircStart(i) <= windsondArr(2).obsTime & windsondArr(2).obsTime <= timeCircEnd(i));
    indWindsond3 = find(timeCircStart(i) <= windsondArr(3).obsTime & windsondArr(3).obsTime <= timeCircEnd(i));
    indWindsond4 = find(timeCircStart(i) <= windsondArr(4).obsTime & windsondArr(4).obsTime <= timeCircEnd(i));
    
    % These arrays represent the times during which the vehicle is flying in
    % a circle.
    procTime1 = windsondArr(1).obsTime(indWindsond1);
    procTime2 = windsondArr(2).obsTime(indWindsond2);
    procTime3 = windsondArr(3).obsTime(indWindsond3);
    procTime4 = windsondArr(4).obsTime(indWindsond4);
    flight(i).times = [procTime1, procTime2, procTime3, procTime4];
    % for plotting temperature vs. time
    flight(i).sensor(1).times = procTime1;
    flight(i).sensor(2).times = procTime2;
    flight(i).sensor(3).times = procTime3;
    flight(i).sensor(4).times = procTime4;
    
    % Circle lats and lons
    latsCirc1 = interp1(iris.obsTimeGPS(indGPS), iris.latitudeGPS_deg(indGPS), procTime1);
    latsCirc2 = interp1(iris.obsTimeGPS(indGPS), iris.latitudeGPS_deg(indGPS), procTime2);
    latsCirc3 = interp1(iris.obsTimeGPS(indGPS), iris.latitudeGPS_deg(indGPS), procTime3);
    latsCirc4 = interp1(iris.obsTimeGPS(indGPS), iris.latitudeGPS_deg(indGPS), procTime4);
    latsCirc = [latsCirc1, latsCirc2, latsCirc3, latsCirc4];

    lonsCirc1 = interp1(iris.obsTimeGPS(indGPS), iris.longitudeGPS_deg(indGPS), procTime1);
    lonsCirc2 = interp1(iris.obsTimeGPS(indGPS), iris.longitudeGPS_deg(indGPS), procTime2);
    lonsCirc3 = interp1(iris.obsTimeGPS(indGPS), iris.longitudeGPS_deg(indGPS), procTime3);
    lonsCirc4 = interp1(iris.obsTimeGPS(indGPS), iris.longitudeGPS_deg(indGPS), procTime4);
    lonsCirc = [lonsCirc1, lonsCirc2, lonsCirc3, lonsCirc4];
    
    % Temps at circle lat/lons
    tempsCirc1 = windsondArr(1).temperature_C(indWindsond1);
    tempsCirc2 = windsondArr(2).temperature_C(indWindsond2);
    tempsCirc3 = windsondArr(3).temperature_C(indWindsond3);
    tempsCirc4 = windsondArr(4).temperature_C(indWindsond4);
    tempsCirc = [tempsCirc1, tempsCirc2, tempsCirc3, tempsCirc4];
    flight(i).temps = tempsCirc;
    % for plotting temperature vs. time
    flight(i).sensor(1).temps = tempsCirc1;
    flight(i).sensor(2).temps = tempsCirc2;
    flight(i).sensor(3).temps = tempsCirc3;
    flight(i).sensor(4).temps = tempsCirc4;
    
    % Find mean temp of each waypoint
    latsMean = [];
    lonsMean = [];
    tempsMean = [];
    for j = 1:nClusters
        indsLat = find(latMin(j) <= latsCirc & latsCirc <= latMax(j));
        indsLon = find(lonMin(j) <= lonsCirc & lonsCirc <= lonMax(j));
        lia = ismember(indsLat,indsLon);
        if isempty(find(lia,1))
            latClust = nan;
            lonClust = nan;
            tempClust = nan;
        else
            loc = find(lia);
            indsClust = indsLat(loc);
            latClust = mean(latsCirc(indsClust));
            lonClust = mean(lonsCirc(indsClust));
            tempClust = mean(tempsCirc(indsClust));
        end
        latsMean = [latsMean latClust];
        lonsMean = [lonsMean lonClust];
        tempsMean = [tempsMean tempClust];
    end
    flight(i).latMean = latsMean;
    flight(i).lonMean = lonsMean;
    flight(i).tempMean = tempsMean;
end
%% Calculate CT^2

% Convert mean lat/lons to mean x/y
for i = 1:nCircles
    [x,y,utmzone] = deg2utm(flight(i).latMean, flight(i).lonMean);
    flight(i).xMean = x;
    flight(i).yMean = y;
end

nRadii = floor(nClusters / 2);
for i = 1:nCircles
    flight(i).radius = struct('length',[],'tempDiff',[]);
    for j = 2:nClusters
        for k = 1:(j-1)
            deltaX = flight(i).xMean(j) - flight(i).xMean(k);
            deltaY = flight(i).yMean(j) - flight(i).yMean(k);
            deltaT = flight(i).tempMean(j) - flight(i).tempMean(k);
            r = sqrt(deltaX^2 + deltaY^2);
            if (j - k) <= nRadii
                n = j - k;
            else
                n = nClusters - (j - k);
            end
            flight(i).radius(n).length(k) = r;
            flight(i).radius(n).tempDiff(k) = deltaT;
        end
    end
    for n = 1:nRadii
        meanLength = mean(flight(i).radius(n).length);
        flight(i).radMean(n) = meanLength;
    end
end


% For plotting temperature difference vs. circumference
baseAngle = 2 * pi / nClusters;
for i = 1:nCircles
    circRad = max(flight(i).radMean) / 2;
    circumference = [];
    for j = 1:nClusters - 1
        circ = circRad * baseAngle * j;
        circumference = [circumference circ];
        tDiffs = [];
        for k = j + 1:nClusters
            dt = sqrt((flight(i).tempMean(k) - flight(i).tempMean(k-j))^2);
            tDiffs = [tDiffs dt];
        end
        flight(i).tempDiffs(j) = mean(tDiffs);
    end
    flight(i).circumference = circumference;
end


DT = [];
for n = 1:nRadii
    sqT = [];
    for i = 1:nCircles
        for j = 1:length(flight(i).radius(n).tempDiff)
            sqTpoint = (flight(i).radius(n).tempDiff(j))^2;
            sqT = [sqT sqTpoint];
        end
    end
    DTpoint = mean(sqT);
    DT = [DT DTpoint];
end

% Calculate CT^2!!!
sepRadius = [];
CT2 = [];
for n = 1:nRadii
    rad = [];
    for i = 1:nCircles
        radPoint = flight(i).radMean(n);
        rad = [rad radPoint];
    end
    sepRadPoint = mean(rad);
    sepRadius = [sepRadius sepRadPoint];
    CT2point = DT(n) / (sepRadPoint^(2/3));
    CT2 = [CT2 CT2point];
end



%% Plots & figures & things

DT_range = [1e-2 1e0];
CT2_range = [1e-3 1e-1];
temp_range = [10 15];

close all

figure(1)
clf
plot(iris.longitudeGPS_deg, iris.latitudeGPS_deg)
xlabel('GPS longitude (deg)')
ylabel('GPS latitude (deg)')
title('IRIS+ 2D Flight Path')
if imgFlag
    imgname = sprintf('./imgs/%4.4d%2.2d%2.2d%s_flightpath2d',...
        procYear,procMonth,procDay,vehicleType);
    print('-f1',imgname,'-dpng');
end

figure(2)
clf
plot3(iris.longitudeGPS_deg, iris.latitudeGPS_deg, iris.altitudeGPS_m)
xlabel('GPS longitude (deg)')
ylabel('GPS latitude (deg)')
zlabel('GPS altitude (m)')
title('IRIS+ 3D Flight Path')
if imgFlag
    imgName = sprintf('./imgs/%4.4d%2.2d%2.2d%s_flightpath3d',...
        procYear,procMonth,procDay,vehicleType);
    print('-f2',imgName,'-dpng');
end

figure(3)
clf
plot(sepRadius, DT)
xlabel('Separation radius (m)')
ylabel('Temperature structure function')
if imgFlag
    imgName = sprintf('./imgs/%4.4d%2.2d%2.2d%s_structurefunction',...
        procYear,procMonth,procDay,vehicleType);
    print('-f3',imgName,'-dpng');
end

figure(4)
clf
semilogy(sepRadius, DT)
set(gca, 'ylim', DT_range)
xlabel('Separation radius (m)')
ylabel('Temperature structure function')
if imgFlag
    imgName = sprintf('./imgs/%4.4d%2.2d%2.2d%s_structurefunctionlog',...
        procYear,procMonth,procDay,vehicleType);
    print('-f4',imgName,'-dpng');
end

figure(5)
clf
plot(sepRadius, CT2)
xlabel('Separation radius (m)')
ylabel('Temperature structure-function parameter')
if imgFlag
    imgName = sprintf('./imgs/%4.4d%2.2d%2.2d%s_structureparam',...
        procYear,procMonth,procDay,vehicleType);
    print('-f5',imgName,'-dpng');
end

figure(6)
clf
semilogy(sepRadius, CT2)
set(gca, 'ylim', CT2_range)
xlabel('Separation radius (m)')
ylabel('Temperature structure-function parameter')
if imgFlag
    imgName = sprintf('./imgs/%4.4d%2.2d%2.2d%s_structureparamlog',...
        procYear,procMonth,procDay,vehicleType);
    print('-f6',imgName,'-dpng');
end

figure(7)
clf
plts = [];
hold on
for n = 1:nSensorFiles
    for i = 1:nCircles
        if ~isnan(windsondArr(n).obsTime(1))
            if n == 1
                color = 'r';
                dispName = sprintf('Windsond %2.2d', n);
            elseif n == 2
                color = 'y';
                dispName = sprintf('Windsond %2.2d', n);
            elseif n == 3
                color = 'g';
                dispName = sprintf('Windsond %2.2d', n);
            elseif n == 4
                color = 'b';
                dispName = sprintf('Windsond %2.2d', n);
            end
            if i == 1
                p = plot(flight(i).sensor(n).times, flight(i).sensor(n).temps, color,...
                    'DisplayName', dispName);
                plts = [plts p];
            else
                plot(flight(i).sensor(n).times, flight(i).sensor(n).temps, color)
            end
            datetick('x',13)
        end
    end
end
hold off
set(gca, 'ylim', temp_range)
xlabel('Observation time')
ylabel('Temperature (C)')
legend(plts,'Location','northwest')
if imgFlag
    imgName = sprintf('./imgs/%4.4d%2.2d%2.2d%s_tempvstime',...
        procYear,procMonth,procDay,vehicleType);
    print('-f7',imgName,'-dpng');
end

figure(8)
clf
hold on
for i = 1:nCircles
    dispName = sprintf('Flight %2.2d', i);
    plot(flight(i).circumference, flight(i).tempDiffs, 'DisplayName', dispName)
end
hold off
xlabel('Circumference (m)')
ylabel('Temperature difference (C)')
legend('show','Location','northwest')
if imgFlag
    imgName = sprintf('./imgs/%4.4d%2.2d%2.2d%s_tempvscircum',...
        procYear,procMonth,procDay,vehicleType);
    print('-f8',imgName,'-dpng');
end
shg




