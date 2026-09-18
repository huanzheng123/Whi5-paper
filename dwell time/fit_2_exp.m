% Plot the histogram and fit the exponential
clear; close all;

% Parameters
interval_time = 1;
truncpt = 2;

% Load data from .mat files
swi4 = load('Swi4_bound_time_results.mat').results;
Whi5 = load('Whi5_bound_time_results.mat').results;
YHZ68 = load('YHZ68_bound_time_results.mat').results;
YHZ75 = load('YHZ75_bound_time_results.mat').results;

% Extract bound time data%
%boundtime_extr1 = [YHZ75{1,1}];
boundtime_extr1 = [YHZ75{4,1}];
%%
% Fit exponential model and estimate bound time
[est, ci, se] = Fitting_truncExponential(boundtime_extr1, interval_time, truncpt);
[Tbound, Tbound_ci, Tbound_err] = Bound_time_estimator_no_bounds(boundtime_extr1, est, 12, 0.1, truncpt);

adjust_fig
axis([0 inf 0 1]);
set(gcf,'position',[200 200 380 320]);

function adjust_fig
set(gca, ...
     'Box'         , 'off'     , ...
     'TickDir'     , 'out'     , ...
     'XMinorTick'  , 'off'      , ...
     'YMinorTick'  , 'off'      , ...
      'YGrid'       , 'off',...
      'XGrid'       , 'off');
set(gca,'LineWidth',1,'FontSize',16,'box','off');
end
