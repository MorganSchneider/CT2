function [ct2vec, ct2av, sepRadius, dAlt, temps, obstime] = CT2_copter(procYear, procMonth, procDay, fileNum, vehicleNum, correction, dAlt, multiple)
% This script will be used to derive the temperature structure-function
% parameter from CopterSonde data.
%
% ***ADD CURVE FIT***
% 
% Credit to Austin Dixon and Phil Chilson for wind estimation code.
% 
% Inputs for debugging:
%
% procYear = 2017;
% procMonth = 12;
% procDay = 19;
% fileNum = 2;
% vehicleNum = 1;
% correction = false;
% dAlt = 15;
% multiple = false;

%% Set up & user inputs

% Flag to decide if image file should be created
imgFlag = false;

% *** You will need to change the baseDir for your computer
% This is where your 'thermo' folder lives
baseDir = './Documents/CLOUDMAP/Code/';

%% Read in the data

% Read the copter data
vehicleType = 'coptersonde';
% WAYPOINT FILES:
% Column 5 = # of seconds at waypoint
% Column 9 = latitude
% Column 10 = longitude
wpFileName = sprintf('%4.4d%2.2d%2.2d_wp.waypoints', procYear, procMonth, procDay);
fp = fopen([ baseDir wpFileName]);
fgetl(fp);
wp = textscan(fp,'%*d %*d %*d %*d %f %*f %*f %*f %f %f %*f %*d');
fclose(fp);
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
sensorFileName = sprintf('%s%d_Data_%4.4d%2.2d%2.2d_%d.csv', vehicleType, vehicleNum,...
    procYear, procMonth, procDay, fileNum);
data = csvread([baseDir sensorFileName], 1, 1);
nSensorFiles = floor(length(data(1,14:end)) / 3);
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
    iMet(i).humidity_percent = data(:,(13 + i)); %columns 14, 15, 16, 17
    iMet(i).RHtemperature_C = data(:,(17 + i)); %columns 18, 19, 20, 21
    iMet(i).temperature_C = data(:,(21 + i)) - 273; %columns 22, 23, 24, 25
end


nNans = 10; % arbitrary
for i = nSensorFiles + 1:4
    iMet(i).temperature_C = nan(1, nNans);
    iMet(i).humidity_percent = nan(1, nNans);
    iMet(i).RHtemperature_C = nan(1, nNans);
end

close all
p = zeros(1,nSensorFiles);
hold on
for n = 1:nSensorFiles
    dn = sprintf('Temp %d', n);
    p(n) = plot(iMet(n).obsTime, iMet(n).temperature_C, 'DisplayName', dn);
end
hold off
xlabel('Obs. Time')
ylabel('Temperature (C)')
legend(p,'Location','best')
shg

sf = 1:nSensorFiles;
x = inputdlg('Enter sensors to ignore, separated by spaces:');
xx = str2num(x{:});

s = sf(ismember(sf,xx) == 0);

close all
clear p

%% Find mean temperature of each waypoint

% Create lat/lon clusters
latWindow = 0.00003; %these are not necessarily universal values
lonWindow = 0.00003;
latMin = desiredLats - latWindow;
lonMin = desiredLons - lonWindow;
latMax = desiredLats + latWindow;
lonMax = desiredLons + lonWindow;
altWindow = 0.5; %meters
altMin = dAlt - altWindow;
altMax = dAlt + altWindow;

% Find temps

% This finds the indices of all times between the circle start and end
% times (basically, the times during which the vehicle is flying in the
% circle of constant height) and puts those indices into their own array.
indFlight = find(altMin <= iMet(1).altitude_m & iMet(1).altitude_m <= altMax);

% These arrays represent the times during which the vehicle is flying in
% a circle.
times = iMet(1).obsTime(indFlight);

% Circle lats and lons
lats = iMet(1).latitude_deg(indFlight);
lons = iMet(1).longitude_deg(indFlight);

% Temps at circle lat/lons
sensor = struct('temps',[]);
for i = 1:length(s)
    sensor(i).temps = iMet(s(i)).temperature_C(indFlight);
end

% Attitude at circle lats/lons
roll = iMet(1).roll_deg(indFlight);
pitch = iMet(1).pitch_deg(indFlight);
yaw = iMet(1).yaw_deg(indFlight);

% Calculate time steps
np = length(times) - 1;
dt = zeros(1, np);
for j = 1:np
    dt(j) = etime(datevec(times(j+1)), datevec(times(j)));
end
dtMean = mean(dt); %seconds

% Ask Brian for his Coptersonde wind estimation code
if correction
    nVals = length(roll);
    psi_deg = zeros(nVals, 1);
    az_deg = zeros(nVals, 1);
    for j = 1:nVals
        crol = cosd(roll(j));
        srol = sind(roll(j));
        cpit = cosd(pitch(j));
        spit = sind(pitch(j));
        cyaw = cosd(yaw(j));
        syaw = sind(yaw(j));
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
    
    windSpeedCoeff_mps = 14*ones(size(yaw));
    % Parameter obtained experimentally to calculate the wind speed
    windSpeedCoeff_mps(yaw <= 50 | yaw >= 300) = 14;
    windSpeed_mps = windSpeedCoeff_mps .* sqrt(tand(psi_deg));
    windDir_deg = az_deg;
    
    ind1 = find(windDir_deg < 0);
    ind2 = find(windDir_deg > 360);
    windDir_deg(ind1) = windDir_deg(ind1) + 360;
    windDir_deg(ind2) = windDir_deg(ind2) - 360;
    
    u = windSpeed_mps .* sind(windDir_deg);
    v = windSpeed_mps .* cosd(windDir_deg);
else
    u = zeros(1, length(indFlight));
    v = zeros(1, length(indFlight));
end

[x, y] = deg2utm(lats, lons);
u_copter = zeros(1, length(dt));
v_copter = zeros(1, length(dt));
for i = 1: length(dt)
    u_copter(i) = (x(i+1) - x(i)) / dt(i);
    v_copter(i) = (y(i+1) - y(i)) / dt(i);
end

cluster = struct('inds', []);

% Find mean temp of each waypoint
latMean = zeros(1, nClusters);
lonMean = zeros(1, nClusters);
tempMean = zeros(1, nClusters);
uMean = zeros(1, nClusters);
vMean = zeros(1, nClusters);
dt = zeros(1, nClusters);
for j = 1:nClusters
    %find indices of each cluster
    indsLat = find(latMin(j) <= lats & lats <= latMax(j));
    indsLon = find(lonMin(j) <= lons & lons <= lonMax(j));
    lia = ismember(indsLat,indsLon);
    %returns array containing logical 1 where data from indsLon is found in indsLat;
    %otherwise, 0
    %lia = "location in A"
    if isempty(find(lia,1))
        % if there is no data within the window, cluster values = nan
        latClust = nan;
        lonClust = nan;
        tempClust = nan;
        uClust = nan;
        vClust = nan;
        dTimes = nan;
    else
        indsClust = indsLat(lia ~= 0);
        % differentiate between first and last waypoint, since they are the
        % exact same lat/lon coordinates
        if (j == 1) || (j == nClusters) %first or last waypoint
            for m = 1:length(indsClust) - 1 %for length of list of indices
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
        cluster(j).inds = indsClust;
        
        % mean values for each cluster
        latClust = nanmean(lats(indsClust));
        lonClust = nanmean(lons(indsClust));
        uClust = nanmean(u(indsClust));
        vClust = nanmean(v(indsClust));
        dTimes = length(indsClust) * dtMean;
        sensorMeanTemp = zeros(1, nSensorFiles);
        %mean temperature at waypoint from each sensor
        for n = 1:length(s)
            sensorMeanTemp(n) = nanmean(sensor(n).temps(indsClust));
        end
        %overall cluster mean temp
        tempClust = nanmean(sensorMeanTemp);
    end
    latMean(j) = latClust;
    lonMean(j) = lonClust;
    tempMean(j) = tempClust;
    uMean(j) = uClust;
    vMean(j) = vClust;
    dt(j) = dTimes;
end

%% Calculate CT^2

% Convert mean lat/lons to mean x/y
[xMean,yMean] = deg2utm(latMean, lonMean);

%advection correction in calculation of separation radii
nRadii = floor(nClusters / 2);
radius = struct('length',[],'tempDiff',[]);
for j = 2:nClusters
    for k = 1:(j-1)
        adx = xMean(j) - xMean(k);
        ady = yMean(j) - yMean(k);
        udt = uMean(k) * dt(k);
        vdt = vMean(k) * dt(k);
        
        %Re = sqrt((x2-x1-udt)^2 + (y2-y1-vdt)^2)
        if (j - k) <= nRadii
            n = j - k;
        else
            n = nClusters - (j - k);
        end
        radius(n).length(k) = sqrt((adx - udt)^2 + (ady - vdt)^2);
        radius(n).tempDiff(k) = tempMean(j) - tempMean(k);
    end
end
radMean = zeros(1, nRadii);
for n = 1:nRadii
    radMean(n) = nanmean(radius(n).length);
end

% For plotting temperature difference vs. circumference
baseAngle = 2 * pi / nClusters;

circRad = max(radMean) / 2;
circumference = zeros(1, nClusters - 1);
tempDiffs = zeros(1, nClusters - 1);
for j = 1:nClusters - 1
    circ = circRad * baseAngle * j;
    circumference(j) = circ;
    tDiffs = zeros(1, nClusters - j);
    %mean temperature variance
    for k = j + 1:nClusters
        tDiffs(k) = sqrt((tempMean(k) - tempMean(k-j))^2);
    end
    tempDiffs(j) = nanmean(tDiffs);
end


% Calculate D(T)
DT = zeros(1, nRadii);
for n = 1:nRadii
    sqT = zeros(1, length(radius(n).tempDiff));
    for j = 1:length(radius(n).tempDiff)
        sqT(j) = (radius(n).tempDiff(j))^2;
    end
    DT(n) = nanmean(sqT);
end

% Calculate CT^2!!!
sepRadius = zeros(1, nRadii);
CT2 = zeros(1, nRadii);
for n = 1:nRadii
    sepRadius(n) = nanmean(radMean(n));
    CT2(n) = DT(n) / (sepRadius(n)^(2/3));
end
%vector of CT2 values per separation radius
ct2vec = CT2;
%Print mean CT2 value for entire flight
ct2av = mean(CT2);
fprintf('CT2 at %d m = %d \n', dAlt, ct2av);

temps = sensor;
obstime = times;

%% Plots & figures & things

close all

if multiple
    if correction
        qual = sprintf('%dm_corr_%d', dAlt, fileNum);
        imgDir = sprintf('imgs/corrected/');
    else
        qual = sprintf('%dm_%d', dAlt, fileNum);
        imgDir = sprintf('imgs/uncorrected/');
    end
else
    if correction
        qual = sprintf('%dm_corr', dAlt);
        imgDir = sprintf('imgs/corrected/');
    else
        qual = sprintf('%dm', dAlt);
        imgDir = sprintf('imgs/uncorrected/');
    end
end

figure(1);
clf
plot(iMet(1).longitude_deg, iMet(1).latitude_deg);
xlabel('GPS longitude (deg)')
ylabel('GPS latitude (deg)')
title(sprintf('CopterSonde 2D flight path at %d m', dAlt))
if imgFlag
    imgName = sprintf('%4.4d%2.2d%2.2d_flightpath2d_%s',...
        procYear,procMonth,procDay,qual);
    print('-f1',[baseDir imgDir imgName],'-dpng');
end

figure(2);
clf
plot3(iMet(1).longitude_deg, iMet(1).latitude_deg, iMet(1).altitude_m);
xlabel('GPS longitude (deg)')
ylabel('GPS latitude (deg)')
zlabel('GPS altitude (m)')
title(sprintf('CopterSonde 3D flight path at %d m', dAlt))
if imgFlag
    imgName = sprintf('%4.4d%2.2d%2.2d_flightpath3d_%s',...
        procYear,procMonth,procDay,qual);
    print('-f2',[baseDir imgDir imgName],'-dpng');
end

figure(3);
clf
plot(sepRadius, DT);
if correction
    xlabel('Adjusted separation radius (m)')
else
    xlabel('Separation radius (m)')
end
ylabel('Temperature structure function')
title(sprintf('D_{T} at %d m', dAlt))
if imgFlag
    imgName = sprintf('%4.4d%2.2d%2.2d_DT_%s',...
        procYear,procMonth,procDay,qual);
    print('-f3',[baseDir imgDir imgName],'-dpng');
end

figure(4);
clf
semilogy(sepRadius, DT);
y = ylim;
ex1 = floor(log10(y(1))) - 1;
ex2 = ceil(log10(y(2))) + 1;
ax = gca;
ax.YLim = [10^ex1 10^ex2];
if correction
    xlabel('Adjusted separation radius (m)')
else
    xlabel('Separation radius (m)')
end
ylabel('Temperature structure function')
title(sprintf('D_{T} at %d m', dAlt))
if imgFlag
    imgName = sprintf('%4.4d%2.2d%2.2d_logDT_%s',...
        procYear,procMonth,procDay,qual);
    print('-f4',[baseDir imgDir imgName],'-dpng');
end

figure(5);
clf
plot(sepRadius, CT2);
if correction
    xlabel('Adjusted separation radius (m)')
else
    xlabel('Separation radius (m)')
end
ylabel('Temperature structure-function parameter')
title(sprintf('C_{T}^{2} at %d m', dAlt))
if imgFlag
    imgName = sprintf('%4.4d%2.2d%2.2d_CT2_%s',...
        procYear,procMonth,procDay,qual);
    print('-f5',[baseDir imgDir imgName],'-dpng');
end

figure(6);
clf
semilogy(sepRadius, CT2);
y = ylim;
ex1 = floor(log10(y(1))) - 1;
ex2 = ceil(log10(y(2))) + 1;
ax = gca;
ax.YLim = [10^ex1 10^ex2];
if correction
    xlabel('Adjusted separation radius (m)')
else
    xlabel('Separation radius (m)')
end
ylabel('Temperature structure-function parameter')
title(sprintf('C_{T}^{2} at %d m', dAlt))
if imgFlag
    imgName = sprintf('%4.4d%2.2d%2.2d_logCT2_%s',...
        procYear,procMonth,procDay,qual);
    print('-f6',[baseDir imgDir imgName],'-dpng');
end

figure(7);
clf
p = zeros(1,length(s));
hold on
for n = 1:length(s)
    if ~isnan(iMet(s(n)).temperature_C(1))
        if n == 1
            color = 'r';
        elseif n == 2
            color = 'y';
        elseif n == 3
            color = 'g';
        elseif n == 4
            color = 'b';
        end
        dispName = sprintf('%s %d', sensorType, n);
        p(n) = plot(times, sensor(n).temps, color,...
            'DisplayName', dispName);
    end
end
meanInds = zeros(1,nClusters);
clusterMeanTemp = zeros(1,nClusters);
for j = 1:nClusters
    if ~isempty(cluster(j).inds)
        vline([times(cluster(j).inds(1)) times(cluster(j).inds(end))], '-k');
        meanInds(j) = floor(nanmean(cluster(j).inds));
        q = zeros(1,length(s));
        for i = 1:length(s)
            q(i) = nanmean(sensor(i).temps(cluster(j).inds));
        end
        clusterMeanTemp(j) = mean(q);
    else
        meanInds(j)= nan;
        clusterMeanTemp(j) = nan;
    end
end
plot(times(meanInds(~isnan(meanInds))), clusterMeanTemp(~isnan(clusterMeanTemp)), '*k');
hold off
xlabel('Observation time')
ylabel('Temperature (C)')
title(sprintf('Observed temperatures at %d m', dAlt))
legend(p,'Location','best')
legend('boxoff')
datetick('x',15)
if imgFlag
    imgName = sprintf('%4.4d%2.2d%2.2d_tempvstime_%s',...
        procYear,procMonth,procDay,qual);
    print('-f7',[baseDir imgDir imgName],'-dpng');
end


figure(8);
clf
plot(circumference, tempDiffs);
xlabel('Circumference (m)')
ylabel('Temperature difference (C)')
title(sprintf('Temperature differences at %d m', dAlt))
if imgFlag
    imgName = sprintf('%4.4d%2.2d%2.2d_tempvscirc_%s',...
        procYear,procMonth,procDay,qual);
    print('-f8',[baseDir imgDir imgName],'-dpng');
end

end



