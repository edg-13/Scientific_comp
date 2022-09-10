warning off; close all; clear all;
% % Change default axes fonts.
set(0,'DefaultAxesFontName', 'Times New Roman');
set(0,'DefaultAxesFontSize', 20);
% % Change default text fonts.
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 20);
%% Numerical solution
% Load the data and the simulations parameters
data0=load('heateq.dat');
sim_params=load('heateq_outparams');
nx = sim_params(1);
% The data can be reshaped to a more convinient format
data=reshape(data0, nx, size(data0,1)/nx, size(data0,2));
% We can visualise the time evolution of temperature by plotting the
% temperature as a color along the grid and at increasing time
figure;
pcolor(data(:,:,1), data(:,:,2), data(:,:,3));
shading flat; colormap(jet(100)); newcolmap=colormap;
xlabel('Time'); ylabel('Position');
%% Compare numerical and analytic solutions
% Sample analytic solution at these locations and at final time tf=t(end)
xa=0:0.01:1.0;
t=data(1,:,1); tf=t(end);
%Construct analytic solution
ya=sin(pi.*xa).*exp(-pi*pi*tf);
% Plot analytic and numerical solutions
figure;
plot(xa, ya,'k-', 'LineWidth', 2); hold on;
% Assign the spatial grid of numerical solution
xn=data(:,1,2); yn=data(:,end,3);
plot(xn, yn, 'ro', 'LineWidth', 2);