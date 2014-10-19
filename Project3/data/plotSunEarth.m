%% Clean up
clear all;
close all;
clc;
format long;


%% Load in data
data = load('SunEarth.txt');


%% Organize data
xS = data(:,2);
yS = data(:,3);
xE = data(:,4);
yE = data(:,5);


%% Plot data
figure(1);
hold on;
plot(xS,yS,'ro');
plot(xE,yE,'b-');
axis equal;
legend('Sun','Earth');