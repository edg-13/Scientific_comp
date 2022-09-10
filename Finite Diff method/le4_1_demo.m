warning off; close all; clear all;
% % Change default axes fonts.
set(0,'DefaultAxesFontName', 'Times New Roman');
set(0,'DefaultAxesFontSize', 20);
% % Change default text fonts.
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 20);
%% Numerical solution
% Load the data and the simulations parameters
g=9.80665;
data=load('le3_1.dat');
inp=load('le3_1_input.txt');
x0=inp(1); y0=inp(2); vx0=inp(3); vy0=inp(4);
tn=data(:,1); xn=data(:,2); yn=data(:,3);
% Analytic solution with no air resistance for comparison
xa=x0:0.01:xn(end);
ya=y0+(vy0/vx0).*(xa-x0)-(0.5*g/(vx0*vx0)).*(xa-x0).*(xa-x0);
%% Numerical solution (black) compared with analytical (red)
figure;
plot(xn, yn, 'k-'); hold on;
plot(xa, ya, 'r.');
xlabel('x'); ylabel('y');