% plot the CI curve
% calculate the confidence interval by puting the 
clear;close all;
% read results and final plot
swi4 = load('Swi4_bound_time_results.mat').results;
Whi5 = load('Whi5_bound_time_results.mat').results;
YHZ68 = load("YHZ68_bound_time_results.mat").results;
YHZ75 = load('YHZ75_bound_time_results.mat').results;

%output the CI results for all the tracks 
datasummary1 = outputcsvdwellbleach(swi4);
datasummary2 = outputcsvdwellbleach(YHZ75);
datasummary3 = outputcsvdwellbleach(YHZ68);
datasummary4 = outputcsvdwellbleach(Whi5);
%%
% generate the data for plot 
model_series = [datasummary1(2,1);datasummary2(2,1);datasummary4(2,1);datasummary3(2,1)];
upper_model_error = [datasummary1(2,2);datasummary2(2,2);datasummary4(2,2);datasummary3(2,2)] % Upper error values
upper_model_error = upper_model_error - model_series; 
lower_model_error = upper_model_error;

b = bar(model_series, 'grouped');
hold on

% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(model_series);

% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end

% Plot the error bars
% Loop through each group and bar to plot individual error bars
for i = 1:nbars
    % Calculate the error for each group for the current bar
    error = [lower_model_error(:,i)'; upper_model_error(:,i)'];
    % Plot the error bars
    errorbar(x(i,:), model_series(:,i), error(1,:), error(2,:), 'k', 'linestyle', 'none','LineWidth',2);
end

adjust_fig

ylabel('Residence time (seconds)');
set(gcf,'position',[200 200 380 320])




%%

model_series = 1-[0.40604 0.46037; 0.47317 0.58238]; 
upper_model_error = 1-[0.4205 0.47938; 0.50727 0.62701];  % Upper error values

upper_model_error = upper_model_error - model_series; 
lower_model_error = upper_model_error;

b = bar(model_series, 'grouped');
hold on

% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(model_series);

% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end

% Plot the error bars
% Loop through each group and bar to plot individual error bars
for i = 1:nbars
    % Calculate the error for each group for the current bar
    error = [lower_model_error(:,i)'; upper_model_error(:,i)'];
    % Plot the error bars
    errorbar(x(i,:), model_series(:,i), error(1,:), error(2,:), 'k', 'linestyle', 'none','LineWidth',2);
end

adjust_fig

ylabel('Survival Probability Decline Rate');
set(gcf,'position',[200 200 380 320])


function adjust_fig
set(gca, ...
     'Box'         , 'off'     , ...
     'TickDir'     , 'out'     , ...
     'XMinorTick'  , 'off'      , ...
     'YMinorTick'  , 'off'      , ...
      'YGrid'       , 'off',...
      'XGrid'       , 'off');
set(gca,'LineWidth',3,'FontSize',16,'box','off');
end

function barplottt(datasummary)
data = datasummary(3:4,1);
errhigh = datasummary(3:4,2) - data;
errlow = data - datasummary(3:4,3);
x = 3:4;
bar(x,data)                
hold on

er = errorbar(x,data,errlow,errhigh);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

end
%%
function datasummary = outputcsvdwellbleach(dataset)
filename = sprintf('%s_fitting_dwell_results.csv','swi4');
datasummary = [];
for celliter = [1,2,3,4,5]
[beta1,CI_beta1] = cal_log_exp_dwell(dataset{celliter,1});
datasummary(celliter,1) = beta1;
datasummary(celliter,2) = CI_beta1(1);
datasummary(celliter,3) = CI_beta1(2);
end
csvwrite(filename,datasummary);
end


function out = getVarName(var)
    out = inputname(1);
end


function [tdwell,CI_beta1_e] = cal_log_exp_dwell(datasets)

[f,x] = ecdf(datasets,'Function','survivor','Alpha',0.01,'Bounds','on');

aaa = find(f == 0);
x(aaa) = [];
f(aaa) = [];

f2 = log(f);
aaa = find(f2 == 0);
x(aaa) = [];
f2(aaa) = [];

aaa = find(x>30);
x(aaa) = [];
f2(aaa) = [];

mdl1 = fitlm(x, f2);

beta1 = mdl1.Coefficients.Estimate(2);
SE_beta1 = mdl1.Coefficients.SE(2);

df1 = length(f2) - 2;
t_crit1 = tinv(0.99, df1);

tduration = -1/beta1;
tdwell = 12*tduration/(12-tduration);

CI_beta1 = [(tdwell -  t_crit1*SE_beta1), (tdwell*SE_beta1)];
CI_beta1_e = (CI_beta1);

end


