% CT2_coptersonde(year, month, day, fileNum, vehicleNum, correction, dAlt, multiple)
%   year: year, 4-digit integer value
%   month: month, 2-digit integer value
%   day: date, 2-digit integer value
%   fileNum: flight number, integer value
%   vehicleNum: which Coptersonde, integer value 1:3 (rest in peace, copter 1...)
%   correction: advection correction, logical value (either 'true'/1 or
%       'false'/0)
%   dAlt: preset altitude for autopilot flight, probably integer value in meters
%
% Returns [ct2vec, ct2av, sepRadius, dAlt]
%   ct2vec: vector of CT2 values per radius, class double
%   ct2av: average CT2 value for the flight, class double
%   sepRadius: vector of separation radii (corresponding to ct2vec values), class double
%   dAlt: flight altitude from input arguments, class integer
%
%

baseDir = './Documents/CLOUDMAP/Code/';

%%  10/03/2017, 4 flights

% consider tossing out sensor 4
[ct2vec_1, ct2av_1, seprad1, dAlt1, temps1, times1] = CT2_copter(2017,10,03,1,1,0,15,1); %uncorrected for advection
[ct2vec_1c, ct2av_1c, ~, ~, ~, ~] = CT2_copter(2017,10,03,1,1,1,15,1); %corrected for advection
[ct2vec_2, ct2av_2, seprad2, dAlt2, temps2, times2] = CT2_copter(2017,10,03,2,1,0,15,1);
[ct2vec_2c, ct2av_2c, ~, ~, ~, ~] = CT2_copter(2017,10,03,2,1,1,15,1);
[ct2vec_3, ct2av_3, seprad3, dAlt3, temps3, times3] = CT2_copter(2017,10,03,3,1,0,30,0);
[ct2vec_3c, ct2av_3c, ~, ~, ~, ~] = CT2_copter(2017,10,03,3,1,1,30,0);
[ct2vec_4, ct2av_4, seprad4, dAlt4, temps4, times4] = CT2_copter(2017,10,03,4,1,0,45,0);
[ct2vec_4c, ct2av_4c, ~, ~, ~, ~] = CT2_copter(2017,10,03,4,1,1,45,0);

close all

% Plot: CT2 vs. separation radius with each flight in a different panel
figure(1);

subplot(4,1,1)
semilogy(seprad1, ct2vec_1, '-b', seprad1, ct2vec_1c, '-r');
title('C_{T}^{2} on 10/03/2017', 'FontSize', 11)
y = ylim;
ex1 = floor(log10(y(1))) - 1;
ex2 = ceil(log10(y(2))) + 1;
ax1 = gca;
ax1.XLim = [0 120];
ax1.YLim = [10^ex1 10^ex2];
set(gca, 'FontSize', 8)
legend('Without correction', 'With correction', 'Location', 'best')
legend('boxoff')
ylabel(sprintf('Flight #1, %d m', dAlt1), 'FontSize', 9)

subplot(4,1,2)
semilogy(seprad2, ct2vec_2, '-b', seprad2, ct2vec_2c, '-r');
y = ylim;
ex1 = floor(log10(y(1))) - 1;
ex2 = ceil(log10(y(2))) + 1;
ax2 = gca;
ax2.XLim = [0 120];
ax2.YLim = [10^ex1 10^ex2];
set(gca, 'FontSize', 8)
ylabel(sprintf('Flight #2, %d m', dAlt2), 'FontSize', 9)

subplot(4,1,3)
semilogy(seprad3, ct2vec_3, '-b', seprad3, ct2vec_3c, '-r');
y = ylim;
ex1 = floor(log10(y(1))) - 1;
ex2 = ceil(log10(y(2))) + 1;
ax3 = gca;
ax3.XLim = [0 120];
ax3.YLim = [10^ex1 10^ex2];
set(gca, 'FontSize', 8)
ylabel(sprintf('Flight #3, %d m', dAlt3), 'FontSize', 9)

subplot(4,1,4)
semilogy(seprad4, ct2vec_4, '-b', seprad4, ct2vec_4c, '-r');
y = ylim;
ex1 = floor(log10(y(1))) - 1;
ex2 = ceil(log10(y(2))) + 1;
ax4 = gca;
ax4.XLim = [0 120];
ax4.YLim = [10^ex1 10^ex2];
set(gca, 'FontSize', 8)
xlabel('Separation radius (m)', 'FontSize', 9)
ylabel(sprintf('Flight #4, %d m', dAlt4), 'FontSize', 9)

imgName = 'imgs/20171003_CT2comparison4p';
print('-f1', [baseDir imgName], '-dpng');

%Plot: CT2 vs separation radius, all flights on same axes
figure(2);
semilogy(seprad1, ct2vec_1, '-r', seprad1, ct2vec_1c, '--r',...
    seprad2, ct2vec_2, '-m', seprad2, ct2vec_2c, '--m',...
    seprad3, ct2vec_3, '-g', seprad3, ct2vec_3c, '--g',...
    seprad4, ct2vec_4, '-b', seprad4, ct2vec_4c, '--b');
y = ylim;
ex1 = floor(log10(y(1))) - 1;
ex2 = ceil(log10(y(2))) + 1;
ax = gca;
ax.XLim = [0 120];
ax.YLim = [10^ex1 10^ex2];
set(gca, 'FontSize', 8)
legend('15m flight #1, uncorrected','15m flight #1, corrected','15m flight #2, uncorrected',...
    '15m flight #2, corrected','30m flight #3, uncorrected','30m flight #3, corrected',...
    '45m flight #4, uncorrected','45m flight #4, corrected','Location','southwest')
legend('boxoff')
xlabel('Separation radius (m)', 'FontSize', 9)
ylabel('C_{T}^{2} (K^{2}m^{-2/3})', 'FontSize', 9)
title('C_{T}^{2} on 10/03/2017', 'FontSize', 11)

imgName = 'imgs/20171003_CT2comparison';
print('-f2', [baseDir imgName], '-dpng');

% Plot: CT2 vs height
figure(3);
semilogy(dAlt1, ct2av_1, 'or', dAlt1, ct2av_1c, '*r',...
    dAlt2, ct2av_2, 'om', dAlt2, ct2av_2c, '*m',...
    dAlt3, ct2av_3, 'og', dAlt3, ct2av_3c, '*g',...
    dAlt4, ct2av_4, 'ob', dAlt4, ct2av_4c, '*b');
y = ylim;
ex1 = floor(log10(y(1))) - 1;
ex2 = ceil(log10(y(2))) + 1;
ax = gca;
ax.XLim = [0 60];
ax.YLim = [10^ex1 10^ex2];
set(gca, 'FontSize', 8)
legend('15m flight #1, uncorrected','15m flight #1, corrected','15m flight #2, uncorrected',...
    '15m flight #2, corrected','30m flight #3, uncorrected','30m flight #3, corrected',...
    '45m flight #4, uncorrected','45m flight #4, corrected','Location','northeast')
xlabel('Height a.g.l. (m)', 'FontSize', 9)
ylabel('C_{T}^{2} (K^{2}m^{-2/3})', 'FontSize', 9)
title('C_{T}^{2} on 10/03/2017', 'FontSize', 11)

imgName = 'imgs/20171003_CT2byHeight';
print('-f3', [baseDir imgName], '-dpng');


%% 10/12/2017, 3 flights

[ct2vec_5, ct2av_5, seprad5, dAlt5, temps5, times5] = CT2_copter(2017,10,12,1,1,0,15,0); %uncorrected
[ct2vec_5c, ct2av_5c, ~, ~, ~, ~]         = CT2_copter(2017,10,12,1,1,1,15,0); %corrected
[ct2vec_6, ct2av_6, seprad6, dAlt6, temps6, times6] = CT2_copter(2017,10,12,2,1,0,30,0);
[ct2vec_6c, ct2av_6c, ~, ~, ~, ~]         = CT2_copter(2017,10,12,2,1,1,30,0);
[ct2vec_7, ct2av_7, seprad7, dAlt7, temps7, times7] = CT2_copter(2017,10,12,3,1,0,50,0);
[ct2vec_7c, ct2av_7c, ~, ~, ~, ~]         = CT2_copter(2017,10,12,3,1,1,50,0);

close all

% Plot: CT2 vs. separation radius with each flight in a different panel
figure(1);

subplot(3,1,1)
semilogy(seprad5, ct2vec_5, '-b', seprad5, ct2vec_5c, '-r');
title('C_{T}^{2} on 10/12/2017', 'FontSize', 11)
y = ylim;
ex1 = floor(log10(y(1))) - 1;
ex2 = ceil(log10(y(2))) + 1;
ax1 = gca;
ax1.XLim = [0 120];
ax1.YLim = [10^ex1 10^ex2];
set(gca, 'FontSize', 8)
legend('Without correction', 'With correction', 'Location', 'best')
legend('boxoff')
ylabel(sprintf('Flight #1, %d m', dAlt5), 'FontSize', 9)

subplot(3,1,2)
semilogy(seprad6, ct2vec_6, '-b', seprad6, ct2vec_6c, '-r');
y = ylim;
ex1 = floor(log10(y(1))) - 1;
ex2 = ceil(log10(y(2))) + 1;
ax2 = gca;
ax2.XLim = [0 120];
ax2.YLim = [10^ex1 10^ex2];
set(gca, 'FontSize', 8)
ylabel(sprintf('Flight #2, %d m', dAlt6), 'FontSize', 9)

subplot(3,1,3)
semilogy(seprad7, ct2vec_7, '-b', seprad7, ct2vec_7c, '-r');
y = ylim;
ex1 = floor(log10(y(1))) - 1;
ex2 = ceil(log10(y(2))) + 1;
ax3 = gca;
ax3.XLim = [0 120];
ax3.YLim = [10^ex1 10^ex2];
set(gca, 'FontSize', 8)
xlabel('Separation radius (m)', 'FontSize', 9)
ylabel(sprintf('Flight #3, %d m', dAlt7), 'FontSize', 9)

imgName = 'imgs/20171012_CT2comparison4p';
print('-f1', [baseDir imgName], '-dpng');

%Plot: CT2 vs separation radius, all flights on same axes
figure(2);
semilogy(seprad5, ct2vec_5, '-r', seprad5, ct2vec_5c, '--r',...
    seprad6, ct2vec_6, '-g', seprad6, ct2vec_6c, '--g',...
    seprad7, ct2vec_7, '-b', seprad7, ct2vec_7c, '--b');
y = ylim;
ex1 = floor(log10(y(1))) - 1;
ex2 = ceil(log10(y(2))) + 1;
ax = gca;
ax.XLim = [0 120];
ax.YLim = [10^ex1 10^ex2];
set(gca, 'FontSize', 8)
legend('15m flight #1, uncorrected','15m flight #1, corrected','30m flight #2, uncorrected',...
    '30m flight #2, corrected','50m flight #3, uncorrected','50m flight #3, corrected',...
    'Location','southwest')
legend('boxoff')
xlabel('Separation radius (m)', 'FontSize', 9)
ylabel('C_{T}^{2} (K^{2}m^{-2/3})', 'FontSize', 9)
title('C_{T}^{2} on 10/12/2017', 'FontSize', 11)

imgName = 'imgs/20171012_CT2comparison';
print('-f2', [baseDir imgName], '-dpng');

% Plot: CT2 vs height
figure(3);
semilogy(dAlt5, ct2av_5, 'or', dAlt5, ct2av_5c, '*r',...
    dAlt6, ct2av_6, 'og', dAlt6, ct2av_6c, '*g',...
    dAlt7, ct2av_7, 'ob', dAlt7, ct2av_7c, '*b');
y = ylim;
ex1 = floor(log10(y(1))) - 1;
ex2 = ceil(log10(y(2))) + 1;
ax = gca;
ax.XLim = [0 60];
ax.YLim = [10^ex1 10^ex2];
set(gca, 'FontSize', 8)
legend('15m flight #1, uncorrected','15m flight #1, corrected','30m flight #2, uncorrected',...
    '30m flight #2, corrected','50m flight #3, uncorrected','50m flight #3, corrected',...
    'Location','northeast')
xlabel('Height a.g.l. (m)', 'FontSize', 9)
ylabel('C_{T}^{2} (K^{2}m^{-2/3})', 'FontSize', 9)
title('C_{T}^{2} on 10/12/2017', 'FontSize', 11)

imgName = 'imgs/20171012_CT2byHeight';
print('-f3', [baseDir imgName], '-dpng');

%% 12/19/2017, 4 flights

% toss out sensor 2
[ct2vec_8, ct2av_8, seprad8, dAlt8, temps8, times8] = CT2_copter(2017,12,19,1,1,0,5,0); %toss out 2 and 4
[ct2vec_8c, ct2av_8c, ~, ~, ~, ~] = CT2_copter(2017,12,19,1,1,1,5,0);
[ct2vec_9, ct2av_9, seprad9, dAlt9, temps9, times9] = CT2_copter(2017,12,19,2,1,0,15,0); %toss out 1 and 2
[ct2vec_9c, ct2av_9c, ~, ~, ~, ~] = CT2_copter(2017,12,19,2,1,1,15,0);
% 30 m data are all incorrect and unusable.
[ct2vec_10, ct2av_10, seprad10, dAlt10, temps10, times10] = CT2_copter(2017,12,19,3,1,0,30,0);
[ct2vec_10c, ct2av_10c, ~, ~, ~, ~] = CT2_copter(2017,12,19,3,1,1,30,0);
[ct2vec_11, ct2av_11, seprad11, dAlt11, temps11, times11] = CT2_copter(2017,12,19,4,1,0,50,0); %toss out 2
[ct2vec_11c, ct2av_11c, ~, ~, ~, ~] = CT2_copter(2017,12,19,4,1,1,50,0);

close all

% Plot: CT2 vs. separation radius with each flight in a different panel
figure(1);

subplot(3,1,1)
semilogy(seprad8, ct2vec_8, '-b', seprad8, ct2vec_8c, '-r');
title('C_{T}^{2} on 12/19/2017', 'FontSize', 11)
y = ylim;
ex1 = floor(log10(y(1))) - 1;
ex2 = ceil(log10(y(2))) + 1;
ax1 = gca;
ax1.XLim = [0 120];
ax1.YLim = [10^ex1 10^ex2];
set(gca, 'FontSize', 8)
legend('Without correction', 'With correction', 'Location', 'best')
legend('boxoff')
ylabel(sprintf('Flight #1, %d m', dAlt8), 'FontSize', 9)

subplot(3,1,2)
semilogy(seprad9, ct2vec_9, '-b', seprad9, ct2vec_9c, '-r');
y = ylim;
ex1 = floor(log10(y(1))) - 1;
ex2 = ceil(log10(y(2))) + 1;
ax2 = gca;
ax2.XLim = [0 120];
ax2.YLim = [10^ex1 10^ex2];
set(gca, 'FontSize', 8)
ylabel(sprintf('Flight #2, %d m', dAlt9), 'FontSize', 9)

subplot(3,1,3)
semilogy(seprad11, ct2vec_11, '-b', seprad11, ct2vec_11c, '-r');
y = ylim;
ex1 = floor(log10(y(1))) - 1;
ex2 = ceil(log10(y(2))) + 1;
ax4 = gca;
ax4.XLim = [0 120];
ax4.YLim = [10^ex1 10^ex2];
set(gca, 'FontSize', 8)
xlabel('Separation radius (m)', 'FontSize', 9)
ylabel(sprintf('Flight #4, %d m', dAlt11), 'FontSize', 9)

imgName = 'imgs/20171219_CT2comparison3p';
%print('-f1', [baseDir imgName], '-dpng');

%Plot: CT2 vs separation radius, all flights on same axes
figure(2);
semilogy(seprad8, ct2vec_8, '-k', seprad8, ct2vec_8c, '--k',...
    seprad9, ct2vec_9, '-r', seprad9, ct2vec_9c, '--r',...
    seprad11, ct2vec_11, '-b', seprad11, ct2vec_11c, '--b');
y = ylim;
ex1 = floor(log10(y(1))) - 1;
ex2 = ceil(log10(y(2))) + 1;
ax = gca;
ax.XLim = [0 120];
ax.YLim = [10^ex1 10^ex2];
set(gca, 'FontSize', 8)
legend('5m flight #1, uncorrected','5m flight #1, corrected','15m flight #2, uncorrected',...
    '15m flight #2, corrected',...
    '50m flight #4, uncorrected','50m flight #4, corrected','Location','southwest')
legend('boxoff')
xlabel('Separation radius (m)', 'FontSize', 9)
ylabel('C_{T}^{2} (K^{2}m^{-2/3})', 'FontSize', 9)
title('C_{T}^{2} on 12/19/2017', 'FontSize', 11)

imgName = 'imgs/20171219_CT2comparison';
%print('-f2', [baseDir imgName], '-dpng');

% Plot: CT2 vs height
figure(3);
semilogy(dAlt8, ct2av_8, 'ok', dAlt8, ct2av_8c, '*k',...
    dAlt9, ct2av_9, 'or', dAlt9, ct2av_9c, '*r',...
    dAlt11, ct2av_11, 'ob', dAlt11, ct2av_11c, '*b');
y = ylim;
ex1 = floor(log10(y(1))) - 1;
ex2 = ceil(log10(y(2))) + 1;
ax = gca;
ax.XLim = [0 60];
ax.YLim = [10^ex1 10^ex2];
set(gca, 'FontSize', 8)
legend('5m flight #1, uncorrected','5m flight #1, corrected','15m flight #2, uncorrected',...
    '15m flight #2, corrected',...
    '50m flight #4, uncorrected','50m flight #4, corrected','Location','northeast')
xlabel('Height a.g.l. (m)', 'FontSize', 9)
ylabel('C_{T}^{2} (K^{2}m^{-2/3})', 'FontSize', 9)
title('C_{T}^{2} on 12/19/2017', 'FontSize', 11)

imgName = 'imgs/20171219_CT2byHeight';
%print('-f3', [baseDir imgName], '-dpng');

%% Combined data
close all

figure(1)
p1 = semilogy(dAlt1, ct2av_1, 'ob');
hold on
p2 = semilogy(dAlt1, ct2av_1c, '*b');
p3 = semilogy(dAlt2, ct2av_2, 'ob'); p4 = semilogy(dAlt2, ct2av_2c, '*b');
p5 = semilogy(dAlt3, ct2av_3, 'og'); p6 = semilogy(dAlt3, ct2av_3c, '*g');
p7 = semilogy(dAlt4, ct2av_4, 'or'); p8 = semilogy(dAlt4, ct2av_4c, '*r');
p9 = semilogy(dAlt5, ct2av_5, 'ob'); p10 = semilogy(dAlt5, ct2av_5c, '*b');
p11 = semilogy(dAlt6, ct2av_6, 'og'); p12 = semilogy(dAlt6, ct2av_6c, '*g');
p13 = semilogy(dAlt7, ct2av_7, 'om'); p14 = semilogy(dAlt7, ct2av_7c, '*m');
p15 = semilogy(dAlt8, ct2av_8, 'ok'); p16 = semilogy(dAlt8, ct2av_8c, '*k');
p17 = semilogy(dAlt9, ct2av_9, 'ob'); p18 = semilogy(dAlt9, ct2av_9c, '*b');
p19 = semilogy(dAlt11, ct2av_11, 'om'); p20 = semilogy(dAlt11, ct2av_11c, '*m');
hold off
y = ylim;
ex1 = floor(log10(y(1))) - 1;
ex2 = ceil(log10(y(2))) + 1;
ax = gca;
ax.XLim = [0 60];
ax.YLim = [10^ex1 10^ex2];
set(gca, 'FontSize', 6)
% gridLegend([p15 p16 p1 p2 p5 p6 p7 p8 p19 p20],5,{'5m uncorrected','5m corrected',...
%     '15m uncorrected','15m corrected','30m uncorrected','30m corrected',...
%     '45m uncorrected','45m corrected','50m uncorrected','50m corrected'},...
%     'Location','southoutside')
xlabel('Height a.g.l. (m)', 'FontSize', 9)
ylabel('C_{T}^{2} (K^{2}m^{-2/3})', 'FontSize', 9)
title('Average C_{T}^{2} on 10/3, 10/12, and 12/19', 'FontSize', 11)

imgName = 'imgs/combined_CT2byHeight';
print('-f1', [baseDir imgName], '-dpng');

figure(2)
semilogy(seprad1, ct2vec_1, '-ob', seprad1, ct2vec_1c, '-*b',...
    seprad2, ct2vec_2, '-ob', seprad2, ct2vec_2c, '-*b',...
    seprad3, ct2vec_3, '-og', seprad3, ct2vec_3c, '-*g',...
    seprad4, ct2vec_4, '-or', seprad4, ct2vec_4c, '-*r',...
    seprad5, ct2vec_5, '-ob', seprad5, ct2vec_5c, '-*b',...
    seprad6, ct2vec_6, '-og', seprad6, ct2vec_6c, '-*g',...
    seprad7, ct2vec_7, '-om', seprad7, ct2vec_7c, '-*m',...
    seprad8, ct2vec_8, '-ok', seprad8, ct2vec_8c, '-*k',...
    seprad9, ct2vec_9, '-ob', seprad9, ct2vec_9c, '-*b',...
    seprad11, ct2vec_11, '-om', seprad11, ct2vec_11c, '-*m');

y = ylim;
ex1 = floor(log10(y(1))) - 1;
ex2 = ceil(log10(y(2))) + 1;
ax = gca;
ax.XLim = [0 120];
ax.YLim = [10^ex1 10^ex2];
xlabel('Separation radius (m)', 'FontSize', 9)
ylabel('C_{T}^{2} (K^{2}m^{-2/3})', 'FontSize', 9)
title('C_{T}^{2} on 10/3, 10/12, and 12/19', 'FontSize', 11)

imgName = 'imgs/combined_CT2comparison';
print('-f2', [baseDir imgName], '-dpng');

figure(3)

p1 = plot(times1, temps1(1).temps+273, '-r', times1, temps1(2).temps+273, '-g',...
    times1, temps1(3).temps+273, '-y', times1, temps1(4).temps+273, '-b');
hold on
p2 = plot(times2, temps2(1).temps+273, '-r', times2, temps2(2).temps+273, '-g',...
    times2, temps2(3).temps+273, '-y', times2, temps2(4).temps+273, '-b');
p3 = plot(times3, temps3(1).temps+273, '-r', times3, temps3(2).temps+273, '-g',...
    times3, temps3(3).temps+273, '-y', times3, temps3(4).temps+273, '-b');
p4 = plot(times4, temps4(1).temps+273, '-r', times4, temps4(2).temps+273, '-g',...
    times4, temps4(3).temps+273, '-y', times4, temps4(4).temps+273, '-b');
hold off
xlabel('Time UTC')
ylabel('Temperature (C)')
datetick('x', 15)
title('Temp. measurements 10/03')
legend(p1,'iMet 1','iMet 2','iMet 3','iMet 4','Location','best')

imgName = 'imgs/20171003_temps';
print('-f3', [baseDir imgName], '-dpng');

figure(4)

p5 = plot(times5, temps5(1).temps+273, '-r', times5, temps5(2).temps+273, '-g',...
    times5, temps5(3).temps+273, '-y', times5, temps5(4).temps+273, '-b');
hold on
p6 = plot(times6, temps6(1).temps+273, '-r', times6, temps6(2).temps+273, '-g',...
    times6, temps6(3).temps+273, '-y', times6, temps6(4).temps+273, '-b');
p7 = plot(times7, temps7(1).temps+273, '-r', times7, temps7(2).temps+273, '-g',...
    times7, temps7(3).temps+273, '-y', times7, temps7(4).temps+273, '-b');
hold off
xlabel('Time UTC')
ylabel('Temperature (C)')
title('Temp. measurements 10/12')
datetick('x', 15)
legend(p5,'iMet 1','iMet 2','iMet 3','iMet 4','Location','best')

imgName = 'imgs/20171012_temps';
print('-f4', [baseDir imgName], '-dpng');

figure(5)

p8 = plot(times8, temps8(1).temps, '-r', times8, temps8(2).temps, '-g',...
    times8, temps8(3).temps, '-y', times8, temps8(4).temps, '-b');
hold on
p9 = plot(times9, temps9(1).temps, '-r', times9, temps9(2).temps, '-g',...
    times9, temps9(3).temps, '-y', times9, temps9(4).temps, '-b');
p10 = plot(times10, temps10(1).temps, '-r', times10, temps10(2).temps, '-g',...
    times10, temps10(3).temps, '-y', times10, temps10(4).temps, '-b');
p11 = plot(times11, temps11(1).temps, '-r', times11, temps11(2).temps, '-g',...
    times11, temps11(3).temps, '-y', times11, temps11(4).temps, '-b');
hold off
xlabel('Time UTC')
ylabel('Temperature (C)')
title('Temp. measurements 12/19')
datetick('x', 15)
legend(p8,'iMet 1','iMet 2','iMet 3','iMet 4','Location','best')

imgName = 'imgs/20171219_temps';
print('-f5', [baseDir imgName], '-dpng');

%% EC comparison 12/19

load('20171219_EC.mat')

semilogy(meanVals.obsTime, meanVals.CT2_C2pm2d3, '-r')
datetick('x', 15)
xlabel('Obs. Time (UTC)')
ylabel('Mean C_{T}^{2} (K^{2}m^{-2/3})')
title('3-D Sonic Anemometer C_{T}^{2}, 12/19/2017')
shg
