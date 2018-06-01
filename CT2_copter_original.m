% This script will be used to derive the temperature structure-function
% parameter from CopterSonde data.
%
% Credit to Phil Chilson for set-up/user inputs and interactive time range option.
% Credit to Austin Dixon and Phil Chilson for wind estimation code.
%
% Ask Petra for CSAT reading/QC code
% Ask Brian for his wind estimation code for Coptersonde
% Get CSAT read routine from Chilson
% Talk to SoM for travel funding
%
% 

%% Set up & user inputs

clear all
clc

% Enter date of flight
procYear = 2017;
procMonth = 10;
procDay = 03;
fileNum = 1;
vehicleNum = 1;

% Flag to indicate whether advection correction should be performed
correction = true;

% Flag to decide if image file should be created
imgFlag = false;

% *** You will need to change the baseDir for your computer
% This is where your 'thermo' folder lives
baseDir = '/Users/User/Documents/CLOUDMAP/Code';

%% Read in the data

% Create the directory of the matlab library and add it to the path
%libDir = [ baseDir 'Github' filesep 'thermo' filesep 'matlab' filesep ];
%addpath(libDir)

% Read the copter data
vehicleType = 'CopterSonde';
% Find the appropriate directory based on instrument type
dataDirName = getDataDir(baseDir, procYear, procMonth, procDay, vehicleType);
% WAYPOINT FILES:
% Column 5 = # of seconds at waypoint
% Column 9 = latitude
% Column 10 = longitude
% Column 11 = altitude
wpFileName = sprintf('%4.4d%2.2d%2.2d_wp.waypoints', procYear, procMonth, procDay);
fp = fopen([ dataDirName wpFileName]);
fgetl(fp);
wp = textscan(fp,'%*d %*d %*d %*d %f %*f %*f %*f %f %f %f %*d');
wpLoiter = wp{1};
wpLats = wp{2};
wpLons = wp{3};
desiredLats = wpLats(wpLoiter~=0);
desiredLons = wpLons(wpLoiter~=0);
if length(desiredLats) ~= length(desiredLons)
    error('Desired lat and lon vectors are not the same length');
end
nClusters = length(desiredLats);

%-------------------------------------------------------------------

% Read the sensor data
sensorType = 'iMet';
% Find the appropriate directory based on instrument type
dataDirName = getDataDir(baseDir, procYear, procMonth, procDay, sensorType);

sensorFileName = sprintf('%s%d_%4.4d-%2.2d-%2.2d_%d.csv', vehicleType, vehicleNum,...
    procYear, procMonth, procDay, fileNum);
data = csvread([dataDirName sensorFileName], 1, 1);
nSensorFiles = floor(length(data(1,14:end)) / 4);
iMet = struct('obsTime',[],'latitude_deg',[],'longitude_deg',[],'altitude_m',...
    [],'pressure_hPa',[],'roll_deg',[],'pitch_deg',[],'yaw_deg',[],'temperature_C',...
    [],'humidity_percent',[],'RHtemperature_C',[],'voltage_V',[]);
for i = 1:nSensorFiles
    % iMet time is given in Epoch time, i.e. the number of the seconds
    % since 01/01/1970 at 00:00:00. To hell with logic and practicality.
    iMet(i).obsTime = data(:,1) / 8.64e4 + datenum(1970,1,1); % 8.64e4= # of seconds in a day
    iMet(i).latitude_deg = data(:,2);
    iMet(i).longitude_deg = data(:,3);
    iMet(i).altitude_m = data(:,4);
    iMet(i).pressure_hPa = data(:,5);
    iMet(i).roll_deg = data(:,6);
    iMet(i).pitch_deg = data(:,7);
    iMet(i).yaw_deg = data(:,8);
    iMet(i).temperatureRH_C = data(:,(14+4*(i - 1))); %columns 14, 18, 22, 26
    iMet(i).humidity_percent = data(:,(15+4*(i - 1))); %columns 15, 19, 23, 27
    iMet(i).temperature_C = data(:,(16+4*(i - 1))); %columns 16, 20, 24, 28
    iMet(i).voltage_V = data(:,(17+4*(i - 1))); %columns 17, 21, 25, 29
end

nNans = 10; % arbitrary
for i = nSensorFiles + 1:4
    iMet(i).temperature_C = nan(1, nNans);
    iMet(i).humidity_percent = nan(1, nNans);
    iMet(i).RHtemperature_C = nan(1, nNans);
    iMet(i).voltage_V = nan(1, nNans);
end

% Remove the matlab library
%rmpath(libDir)

%% Find time range. (Interactive)

timeBeg = iMet(1).obsTime(1);
timeEnd = iMet(1).obsTime(end);
indGPS = find(iMet(1).obsTime);

% 05/05/2017: 1 flight circle
nCircles = input('Input number of flights to process: ');
figure(1)
clf
plot(iMet(1).obsTime(indGPS), iMet(1).altitude_m(indGPS))
xlabel('Time (UTC)')
ylabel('CopterSonde GPS Altitude (m)')
datetick('x', 15)
fprintf('Click circle start and end times\n')
shg

timeCircStart = zeros(1, nCircles);
timeCircEnd = zeros(1, nCircles);
for i = 1:nCircles
    [x, ~] = ginput(2);
    starts = x(1);
    ends = x(2);
    if (starts < timeBeg || starts > timeEnd), starts = nan; end
    if (ends < timeBeg || ends > timeEnd), ends = nan; end
    timeCircStart(i) = starts;
    timeCircEnd(i) = ends;
end

%% Find mean temperature of each waypoint

% Create lat/lon clusters
latWindow = 0.00003; %these are not necessarily universal values
lonWindow = 0.00003;
latMin = desiredLats - latWindow;
lonMin = desiredLons - lonWindow;
latMax = desiredLats + latWindow;
lonMax = desiredLons + lonWindow;
altWindow = 1.0; %meters
altMin = desiredAlts - altWindow;
altMax = desiredAlts + altWindow;

% Find temps (each iteration finds cluster temps for one flight).
flight = struct('times',[],'lats',[],'lons',[],'latMean',[],'lonMean',[],...
    'tempMean',[],'xMean',[],'yMean',[],'radMean',[],'circumference',[],...
    'tempDiffs',[],'u',[],'v',[],'uMean',[],'vMean',[],'dt',[],'roll',[],...
    'pitch',[],'yaw',[],'wpAlt',[]);

for i = 1:nCircles %for i = 1:nAlts
    % This finds the indices of all times between the circle start and end
    % times (basically, the times during which the vehicle is flying in the
    % circle of constant height) and puts those indices into their own array.
    indFlights = find(altMin(i) <= iMet(1).altitude_m & iMet(1).altitude_m <= altMax(i));
    
    % These arrays represent the times during which the vehicle is flying in
    % a circle.
    flight(i).times = iMet(1).obsTime(indFlights);
    
    % Circle lats and lons
    flight(i).lats = iMet(1).latitude_deg(indFlights);
    flight(i).lons = iMet(1).longitude_deg(indFlights);
    
    % Programmed flight altitude AGL
    flight(i).wpAlt = desiredAlts(i);
    
    % Temps at circle lat/lons
    flight(i).sensor(1).temps = iMet(1).temperature_C(indFlights);
    flight(i).sensor(2).temps = iMet(2).temperature_C(indFlights);
    flight(i).sensor(3).temps = iMet(3).temperature_C(indFlights);
    flight(i).sensor(4).temps = iMet(4).temperature_C(indFlights);
    
    % Attitude at circle lats/lons
    flight(i).roll = iMet(1).roll_deg(indFlights);
    flight(i).pitch = iMet(1).pitch_deg(indFlights);
    flight(i).yaw = iMet(1).yaw_deg(indFlights);
    
    % Calculate time steps
    np = length(flight(i).times) - 1;
    dt = zeros(1, np);
    for j = 1:np
        dt(j) = etime(datevec(flight(i).times(j+1)), datevec(flight(i).times(j)));
    end
    dtMean = mean(dt); %seconds
    
    % Ask Brian for his Coptersonde wind estimation code
    if correction
        nVals = length(flight(i).roll);
        psi_deg = zeros(nVals, 1);
        az_deg = zeros(nVals, 1);
        for j = 1:nVals
            crol = cosd(flight(i).roll(j));
            srol = sind(flight(i).roll(j));
            cpit = cosd(flight(i).pitch(j));
            spit = sind(flight(i).pitch(j));
            cyaw = cosd(flight(i).yaw(j));
            syaw = sind(flight(i).yaw(j));
            Rx = [[1 0 0]; ...
                [0 crol srol]; ...
                [0 -srol crol]];
            Ry = [[cpit 0 -spit]; ...
                [0 1 0]; ...
                [spit 0 cpit]];
            Rz = [[cyaw -syaw 0]; ...
                [syaw cyaw 0]; ...
                [0 0 1]];
            R = Rz*Ry*Rx;
            vectorRot = R*[0; 0; 1];
            % Inclination angle
            psi_deg(j) = acosd(dot([0; 0; 1], vectorRot));
            az_deg(j) = atan2d(vectorRot(2), vectorRot(1));
        end
        if az_deg < 0, az_deg = az_deg + 360; end
        
        windSpeedCoeff_mps = 14*ones(size(flight(i).yaw));
        % Parameter obtained experimentally to calculate the wind speed
        ind = find(flight(i).yaw <= 50 | flight(i).yaw >= 300);
        windSpeedCoeff_mps(ind) = 14;
        windSpeed_mps = windSpeedCoeff_mps .* sqrt(tand(psi_deg));
        windDir_deg = az_deg;
        
        ind1 = find(windDir_deg < 0);
        ind2 = find(windDir_deg > 360);
        windDir_deg(ind1) = windDir_deg(ind1) + 360;
        windDir_deg(ind2) = windDir_deg(ind2) - 360;
        
        flight(i).u = windSpeed_mps .* sind(windDir_deg);
        flight(i).v = windSpeed_mps .* cosd(windDir_deg);
    else
        flight(i).u = zeros(1, length(indFlights));
        flight(i).v = zeros(1, length(indFlights));
    end
    
    % Find mean temp of each waypoint
    latsMean = zeros(1, nClusters);
    lonsMean = zeros(1, nClusters);
    tempsMean = zeros(1, nClusters);
    uMean = zeros(1, nClusters - 1);
    vMean = zeros(1, nClusters - 1);
    dTime = zeros(1, nClusters - 1);
    for j = 1:nClusters
        indsLat = find(latMin(j) <= flight(i).lats & flight(i).lats <= latMax(j));
        indsLon = find(lonMin(j) <= flight(i).lons & flight(i).lons <= lonMax(j));
        lia = ismember(indsLat,indsLon);
        %returns array containing logical 1 where data in A is found in B;
        %otherwise, 0
        if isempty(find(lia,1))
            latClust = nan;
            lonClust = nan;
            tempClust = nan;
            uClust = nan;
            vClust = nan;
        else
            if j == 1
                f = 1;
            else
                f = indsClust(end) + 1;
            end
            loc = find(lia);
            indsClust = indsLat(loc);
            if (j == 1) || (j == nClusters) %first or last waypoint
                if (desiredLats(1) == desiredLats(end)) && (desiredLons(1) == desiredLons(end))
                    %if first and last waypoints are the same
                    for m = 1:length(indsClust) %for length of list of indices
                        if indsClust(m+1) ~= indsClust(m) + 1 %finds time "split" in indices
                            if j == 1 %eliminate last waypoint data from first waypoint indices
                                indsClust(m+1:end) = [];
                            elseif j == nClusters %eliminate first waypoint data from last waypoint indices
                                indsClust(1:m) = [];
                            end
                            break
                        end
                    end
                end
            end
            l = indsClust(end);
            
            latClust = mean(flight(i).lats(indsClust));
            lonClust = mean(flight(i).lons(indsClust));
            uClust = mean(flight(i).u(f:l)); %need to change this to account for
            vClust = mean(flight(i).v(f:l)); %advection between waypoints, not just at waypoints
            dTimes = length(f:l) * dtMean;
            sensorMeanTemps = zeros(1, nSensorFiles);
            for n = 1:nSensorFiles
                sensorTemps = zeros(1, length(indsClust));
                for m = 1:length(indsClust)
                    if ~isnan(flight(i).sensor(n).temps(m))
                        sensorTemps(m) = flight(i).sensor(n).temps(m);
                    end
                end
                sensorMeanTemps(n) = nanmean(sensorTemps);
            end
            tempClust = mean(sensorMeanTemps);
        end
        latsMean(j) = latClust;
        lonsMean(j) = lonClust;
        tempsMean(j) = tempClust;
        uMean(j) = uClust;
        vMean(j) = vClust;
        dTime(j) = dTimes;
    end
    flight(i).latMean = latsMean;
    flight(i).lonMean = lonsMean;
    flight(i).tempMean = tempsMean;
    flight(i).uMean = uMean;
    flight(i).vMean = vMean;
    flight(i).dt = dTime;
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
            adx = flight(i).xMean(j) - flight(i).xMean(k);
            ady = flight(i).yMean(j) - flight(i).yMean(k);
            deltaT = flight(i).tempMean(j) - flight(i).tempMean(k);
            udt = flight(i).uMean(k) * flight(i).dt(k);
            vdt = flight(i).vMean(k) * flight(i).dt(k);
            %Re = sqrt((x2-x1-udt)^2 + (y2-y1-vdt)^2)
            r = sqrt((adx - udt)^2 + (ady - vdt)^2);
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
    circumference = zeros(1, nClusters - 1);
    for j = 1:nClusters - 1
        circ = circRad * baseAngle * j;
        circumference(j) = circ;
        tDiffs = zeros(1, nClusters - j);
        for k = j + 1:nClusters
            dT = sqrt((flight(i).tempMean(k) - flight(i).tempMean(k-j))^2);
            tDiffs(k) = dT;
        end
        flight(i).tempDiffs(j) = mean(tDiffs);
    end
    flight(i).circumference = circumference;
end


DT = zeros(1, nRadii);
for n = 1:nRadii
    sqT = [];
    for i = 1:nCircles
        for j = 1:length(flight(i).radius(n).tempDiff)
            sqTpoint = (flight(i).radius(n).tempDiff(j))^2;
            sqT = [sqT sqTpoint];
        end
    end
    DTpoint = mean(sqT);
    DT(n) = DTpoint;
end

% Calculate CT^2!!!
sepRadius = zeros(1, nRadii);
CT2 = zeros(1, nRadii);
for n = 1:nRadii
    rad = zeros(1, nCircles);
    for i = 1:nCircles
        radPoint = flight(i).radMean(n);
        rad(i) = radPoint;
    end
    sepRadPoint = mean(rad);
    sepRadius(n) = sepRadPoint;
    CT2point = DT(n) / (sepRadPoint^(2/3));
    CT2(n) = CT2point;
end
CT2_mean = mean(CT2);
fprintf('CT2 = %d', CT2_mean);

%% Plots & figures & things

DT_range = [1e-4 1e-2];
CT2_range = [1e-5 1e-3];

close all

figure(1)
clf
plot(iMet(1).longitude_deg, iMet(1).latitude_deg)
xlabel('GPS longitude (deg)')
ylabel('GPS latitude (deg)')
title('CopterSonde 2D Flight Path')
if imgFlag
    imgname = sprintf('./imgs/%4.4d%2.2d%2.2d%s_flightpath2d',...
        procYear,procMonth,procDay,vehicleType);
    print('-f1',imgname,'-dpng');
end

figure(2)
clf
plot3(iMet(1).longitude_deg, iMet(1).latitude_deg, iMet(1).altitude_m)
xlabel('GPS longitude (deg)')
ylabel('GPS latitude (deg)')
zlabel('GPS altitude (m)')
title('CopterSonde 3D Flight Path')
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
plts = zeros(1, nSensorFiles);
hold on
for n = 1:nSensorFiles
    for i = 1:nCircles
        if ~isnan(iMet(n).obsTime(1))
            if n == 1
                color = 'r';
                dispName = sprintf('%s %2.2d', sensorType, n);
            elseif n == 2
                color = 'y';
                dispName = sprintf('%s %2.2d', sensorType, n);
            elseif n == 3
                color = 'g';
                dispName = sprintf('%s %2.2d', sensorType, n);
            elseif n == 4
                color = 'b';
                dispName = sprintf('%s %2.2d', sensorType, n);
            end
            if i == 1
                p = plot(flight(i).times, flight(i).sensor(n).temps, color,...
                    'DisplayName', dispName);
                plts(i) = p;
            else
                plot(flight(i).times, flight(i).sensor(n).temps, color)
            end
        end
    end
end
hold off
xlabel('Observation time')
ylabel('Temperature (C)')
legend(plts,'Location','northwest')
datetick('x',15)
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
if nCircles > 1
    legend('show','Location','northwest')
end
if imgFlag
    imgName = sprintf('./imgs/%4.4d%2.2d%2.2d%s_tempvscircum',...
        procYear,procMonth,procDay,vehicleType);
    print('-f8',imgName,'-dpng');
end

shg




