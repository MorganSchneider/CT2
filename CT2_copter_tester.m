% This script will be used to derive the temperature structure-function
% parameter from CopterSonde data.
% 
% Credit to Phil Chilson for set-up/user inputs and interactive time range option.
% Credit to Austin Dixon and Phil Chilson for wind estimation code.
% 
% Ask Brian for his wind estimation code for Coptersonde
% Jeff Basara for KAEFS flux station data processing (QC, etc)
%
% Sensor 1: RH and T are switched???
% All sensor structures contain the same data??
% 
% 
% 
% 

%% Set up & user inputs

clear all
clc

% Enter date of flight
procYear = 2017;
procMonth = 10;
procDay = 03;
vehicleNum = 1;
alts = [15 15 30 45];


% Flag to indicate whether advection correction should be performed
correction = true;

% Flag to decide if image file should be created
imgFlag = false;

% Just because I haven't fixed my plotting code yet
plotFlag = false;

% *** You will need to change the baseDir for your computer
% This is where your 'thermo' folder lives
baseDir = '/Users/Morgan/Documents/CLOUDMAP/Code/';

%% Read in the data

% Read the copter data
vehicleType = 'Coptersonde';
% WAYPOINT FILES:
% Column 9 = latitude
% Column 10 = longitude
% Column 11 = altitude
wpFileName = sprintf('%4.4d%2.2d%2.2d_wp.waypoints', procYear, procMonth, procDay);
fp = fopen([ baseDir wpFileName]);
fgetl(fp);
wp = textscan(fp,'%*d %*d %*d %*d %*f %*f %*f %*f %f %f %f %*d');
wpLats_cell = wp{1};
wpLons_cell = wp{2};
wpAlts_cell = wp{3};

wpLats = wpLats_cell(wpAlts_cell ~= 0);
wpLons = wpLons_cell(wpAlts_cell ~= 0);

lats = factorVec(wpLats);
lons = factorVec(wpLons);

desiredAlts = factorVec(alts);
desiredLats = [lats lats(1)];
desiredLons = [lons lons(1)];
if length(desiredLats) ~= length(desiredLons)
    error('Desired lat and lon vectors are not the same length');
end
nAlts = length(desiredAlts);
nClusters = length(desiredLats);

%-------------------------------------------------------------------

% Read the sensor data
sensorType = 'iMet';
% Find the appropriate directory based on instrument type
dataDirName = getDataDir(baseDir, procYear, procMonth, procDay, sensorType);
iMet = struct('obsTime',[],'latitude_deg',[],'longitude_deg',[],'altitude_m',...
    [],'pressure_hPa',[],'roll_deg',[],'pitch_deg',[],'yaw_deg',[],'temperature_C',...
    [],'humidity_percent',[],'RHtemperature_C',[],'voltage_V',[]);
nFlights = 1:length(alts);
for m = 1:length(alts)
    sensorFileName = sprintf('%s%d_Data_%4.4d-%2.2d-%2.2d_%d.csv', vehicleType, vehicleNum,...
        procYear, procMonth, procDay, nFlights(m));
    data = csvread([baseDir sensorFileName], 1, 1);
    nSensorFiles = floor(length(data(1,14:end)) / 4);
    for i = 1:nSensorFiles
        % iMet time is given in Epoch time, i.e. the number of the seconds
        % since 01/01/1970 at 00:00:00. To hell with logic and practicality.
        flight(m).iMet(i).obsTime = rot90(data(:,1) / 8.64e4 + datenum(1970,1,1)); % 8.64e4= # of seconds in a day
        flight(m).iMet(i).latitude_deg = rot90(data(:,2));
        flight(m).iMet(i).longitude_deg = rot90(data(:,3));
        flight(m).iMet(i).altitude_m = rot90(data(:,4));
        %iMet(i).pressure_hPa = data(:,5);
        flight(m).iMet(i).roll_deg = rot90(data(:,6));
        flight(m).iMet(i).pitch_deg = rot90(data(:,7));
        flight(m).iMet(i).yaw_deg = rot90(data(:,8));
        %iMet(i).temperatureRH_C = data(:,(14+4*(i - 1))); %columns 14, 18, 22, 26
        %iMet(i).humidity_percent = data(:,(15+4*(i - 1))); %columns 15, 19, 23, 27
        flight(m).iMet(i).temperature_C = rot90(data(:,(16+4*(i - 1)))); %columns 16, 20, 24, 28
        %iMet(i).voltage_V = data(:,(17+4*(i - 1))); %columns 17, 21, 25, 29
    end


    nNans = 10; % arbitrary
    for i = nSensorFiles + 1:4
        flight(m).iMet(i).temperature_C = nan(1, nNans);
        %iMet(i).humidity_percent = nan(1, nNans);
        %iMet(i).RHtemperature_C = nan(1, nNans);
        %iMet(i).voltage_V = nan(1, nNans);
    end
end

for i = 1:nSensorFiles
    iMetTimes = [];
    iMetLats = [];
    iMetLons = [];
    iMetAlts = [];
    iMetRoll = [];
    iMetPitch = [];
    iMetYaw = [];
    iMetTemps = [];
    for m = 1:length(nFlights)
        iMetTimes = [iMetTimes flight(m).iMet(i).obsTime];
        iMetLats = [iMetLats flight(m).iMet(i).latitude_deg];
        iMetLons = [iMetLons flight(m).iMet(i).longitude_deg];
        iMetAlts = [iMetAlts flight(m).iMet(i).altitude_m];
        iMetRoll = [iMetRoll flight(m).iMet(i).roll_deg];
        iMetPitch = [iMetPitch flight(m).iMet(i).pitch_deg];
        iMetYaw = [iMetYaw flight(m).iMet(i).yaw_deg];
        iMetTemps = [iMetTemps flight(m).iMet(i).temperature_C];
    end
    iMet(i).obsTime = iMetTimes;
    iMet(i).latitude_deg = iMetLats;
    iMet(i).longitude_deg = iMetLons;
    iMet(i).altitude_m = iMetAlts;
    iMet(i).roll_deg = iMetRoll;
    iMet(i).pitch_deg = iMetPitch;
    iMet(i).yaw_deg = iMetYaw;
    iMet(i).temperature_C = iMetTemps;
end

%% Find mean temperature of each waypoint

% Create lat/lon clusters
latWindow = 0.00003; %these are not necessarily universal values
lonWindow = 0.00003;
latMin = desiredLats - latWindow;
lonMin = desiredLons - lonWindow;
latMax = desiredLats + latWindow;
lonMax = desiredLons + lonWindow;
altWindow = 0.5; %meters
altMin = desiredAlts - altWindow;
altMax = desiredAlts + altWindow;

% Find temps (each iteration finds cluster temps for one flight).
height = struct('alt',[],'times',[],'lats',[],'lons',[],'latMean',[],'lonMean',[],...
    'tempMean',[],'xMean',[],'yMean',[],'radMean',[],'circumference',[],...
    'tempDiffs',[],'u',[],'v',[],'uMean',[],'vMean',[],'dt',[],'roll',[],...
    'pitch',[],'yaw',[],'DT',[],'CT2',[]);

for i = 1:nAlts
    % This finds the indices of all times between the circle start and end
    % times (basically, the times during which the vehicle is flying in the
    % circle of constant height) and puts those indices into their own array.
    indFlights = [];
    for k = 1:nSensorFiles
        inds = find(altMin(i) <= iMet(k).altitude_m & iMet(k).altitude_m <= altMax(i));
        indFlights = [indFlights inds];
    end

    
    % These arrays represent the times during which the vehicle is flying in
    % a circle.
    height(i).times = iMet(1).obsTime(indFlights);
    
    % Circle lats and lons
    height(i).lats = iMet(1).latitude_deg(indFlights);
    height(i).lons = iMet(1).longitude_deg(indFlights);
    
    % Programmed flight altitude AGL
    height(i).alt = desiredAlts(i);
    
    % Temps at circle lat/lons
    for k = 1:nSensorFiles
        height(i).sensor(k).temps = iMet(k).temperature_C(indFlights);
    end
    
    % Attitude at circle lats/lons
    height(i).roll = iMet(1).roll_deg(indFlights);
    height(i).pitch = iMet(1).pitch_deg(indFlights);
    height(i).yaw = iMet(1).yaw_deg(indFlights);
    
    % Calculate time steps
    np = length(height(i).times) - 1;
    dt = zeros(1, np);
    for j = 1:np
        dt(j) = etime(datevec(height(i).times(j+1)), datevec(height(i).times(j)));
        dt(dt > 30) = nan;
    end
    dtMean = nanmean(dt); %seconds
    
    % Ask Brian for his Coptersonde wind estimation code
    if correction
        nVals = length(height(i).roll);
        psi_deg = zeros(nVals, 1);
        az_deg = zeros(nVals, 1);
        for j = 1:nVals
            crol = cosd(height(i).roll(j));
            srol = sind(height(i).roll(j));
            cpit = cosd(height(i).pitch(j));
            spit = sind(height(i).pitch(j));
            cyaw = cosd(height(i).yaw(j));
            syaw = sind(height(i).yaw(j));
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
        
        windSpeedCoeff_mps = 14*ones(size(height(i).yaw));
        % Parameter obtained experimentally to calculate the wind speed
        ind = find(height(i).yaw <= 50 | height(i).yaw >= 300);
        windSpeedCoeff_mps(ind) = 14;
        windSpeed_mps = windSpeedCoeff_mps * sqrt(tand(psi_deg));
        windDir_deg = az_deg;
        
        ind1 = find(windDir_deg < 0);
        ind2 = find(windDir_deg > 360);
        windDir_deg(ind1) = windDir_deg(ind1) + 360;
        windDir_deg(ind2) = windDir_deg(ind2) - 360;
        
        height(i).u = windSpeed_mps .* sind(windDir_deg);
        height(i).v = windSpeed_mps .* cosd(windDir_deg);
    else
        height(i).u = zeros(1, length(indFlights));
        height(i).v = zeros(1, length(indFlights));
    end
    
    % Find mean temp of each waypoint
    latsMean = zeros(1, nClusters);
    lonsMean = zeros(1, nClusters);
    tempsMean = zeros(1, nClusters);
    uMean = zeros(1, nClusters - 1);
    vMean = zeros(1, nClusters - 1);
    dTime = zeros(1, nClusters - 1);
    for j = 1:nClusters
        indsLat = find(latMin(j) <= height(i).lats & height(i).lats <= latMax(j));
        indsLon = find(lonMin(j) <= height(i).lons & height(i).lons <= lonMax(j));
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
            loc = find(lia);
            indsClust = indsLat(loc);
            if j == 1
                f = 1;
            else
                f = indsClust(end) + 1;
            end
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
            
            latClust = mean(height(i).lats(indsClust));
            lonClust = mean(height(i).lons(indsClust));
            uClust = mean(height(i).u(f:l));
            vClust = mean(height(i).v(f:l));
            sensorMeanTemps = zeros(1, nSensorFiles);
            for n = 1:nSensorFiles
                sensorTemps = zeros(1, length(indsClust));
                for m = 1:length(indsClust)
                    if ~isnan(height(i).sensor(n).temps(m))
                        sensorTemps(m) = height(i).sensor(n).temps(m);
                    end
                end
                sensorMean = mean(sensorTemps);
                sensorMeanTemps(n) = sensorMean;
            end
            tempClust = mean(sensorMeanTemps);
            dTimes = length(f:l) * dtMean;
            dTime(j) = dTimes;
        end
        latsMean(j) = latClust;
        lonsMean(j) = lonClust;
        tempsMean(j) = tempClust;
        uMean(j) = uClust;
        vMean(j) = vClust;
    end
    height(i).latMean = latsMean;
    height(i).lonMean = lonsMean;
    height(i).tempMean = tempsMean;
    height(i).uMean = uMean;
    height(i).vMean = vMean;
    height(i).dt = dTime;
end

%% Calculate CT^2

% Convert mean lat/lons to mean x/y
for i = 1:nAlts
    [x,y,utmzone] = deg2utm(height(i).latMean, height(i).lonMean);
    height(i).xMean = x;
    height(i).yMean = y;
end

nRadii = floor(nClusters / 2);
for i = 1:nAlts
    height(i).radius = struct('length',[],'tempDiff',[]);
    for j = 2:nClusters
        for k = 1:(j-1)
            adx = height(i).xMean(j) - height(i).xMean(k);
            ady = height(i).yMean(j) - height(i).yMean(k);
            deltaT = height(i).tempMean(j) - height(i).tempMean(k);
            udt = height(i).uMean(k) * height(i).dt(k);
            vdt = height(i).vMean(k) * height(i).dt(k);
            %Re = sqrt((x2-x1-udt)^2 + (y2-y1-vdt)^2)
            r = sqrt((adx - udt)^2 + (ady - vdt)^2);
            if (j - k) <= nRadii
                n = j - k;
            else
                n = nClusters - (j - k);
            end
            height(i).radius(n).length(k) = r;
            height(i).radius(n).tempDiff(k) = deltaT;
        end
    end
    for n = 1:nRadii
        height(i).radmean(n) = mean(height(i).radius(n).length);
    end
end


% For plotting temperature difference vs. circumference
baseAngle = 2 * pi / nClusters;
for i = 1:nAlts
    circRad = max(height(i).radMean) / 2;
    circumference = zeros(1, nClusters - 1);
    for j = 1:nClusters - 1
        circ = circRad * baseAngle * j;
        circumference(j) = circ;
        tDiffs = zeros(1, nClusters - j);
        for k = j + 1:nClusters
            dT = sqrt((height(i).tempMean(k) - height(i).tempMean(k-j))^2);
            tDiffs(k) = dT;
        end
        height(i).tempDiffs(j) = mean(tDiffs);
    end
    height(i).circumference = circumference;
end

for i = 1:nAlts
    height(i).DT = zeros(1, nRadii);
    for n = 1:nRadii
        sqT = [];
        for j = 1:length(height(i).radius(n).tempDiff)
            sqTpoint = (height(i).radius(n).tempDiff(j))^2;
            sqT = [sqT sqTpoint];
        end
        height(i).DT(n) = mean(sqT);
    end
end

% Calculate CT^2!!!
for i = 1:nAlts
    sepRadius = zeros(1, nRadii);
    CT2 = zeros(1, nRadii);
    for n = 1:nRadii
        CT2(n) = DT(n) / (height(i).radmean(n)^(2/3));
    end
    height(i).CT2 = mean(CT2);
    fprintf('height(%d).CT2 = %d', i, CT2_mean);
end

%% Plots & figures & things

DT_range = [1e-4 1e-2];
CT2_range = [1e-5 1e-3];

if plotFlag
close all

figure(1)
clf
plotFlag(iMet(1).longitude_deg, iMet(1).latitude_deg)
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
% PLOT height(i).radmean vs height(i).DT
clf
plotFlag(sepRadius, DT)
xlabel('Separation radius (m)')
ylabel('Temperature structure function')
if imgFlag
    imgName = sprintf('./imgs/%4.4d%2.2d%2.2d%s_structurefunction',...
        procYear,procMonth,procDay,vehicleType);
    print('-f3',imgName,'-dpng');
end

figure(4)
% PLOT height(i).radmean vs height(i).DT
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
% height(i).radmean vs height(i).CT2
clf
plotFlag(sepRadius, CT2)
xlabel('Separation radius (m)')
ylabel('Temperature structure-function parameter')
if imgFlag
    imgName = sprintf('./imgs/%4.4d%2.2d%2.2d%s_structureparam',...
        procYear,procMonth,procDay,vehicleType);
    print('-f5',imgName,'-dpng');
end

figure(6)
% height(i).radmean vs height(i).CT2
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
% this is fine
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
                p = plotFlag(height(i).times, height(i).sensor(n).temps, color,...
                    'DisplayName', dispName);
                plts(i) = p;
            else
                plotFlag(height(i).times, height(i).sensor(n).temps, color)
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
% this is probably also fine
clf
hold on
for i = 1:nCircles
    dispName = sprintf('Flight %2.2d', i);
    plotFlag(height(i).circumference, height(i).tempDiffs, 'DisplayName', dispName)
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

end



