% calculate the confidence interval by puting the 

% read results and final plot
swi4 = load('Swi4_bound_time_results.mat').results;
Whi5 = load('Whi5_bound_time_results.mat').results;
YHZ68 = load("YHZ68_bound_time_results.mat").results;
YHZ75 = load('YHZ75_bound_time_results.mat').results;

% output the CI results for all the tracks 
% outputcsvdwell(swi4)
% outputcsvdwell(Whi5)
% outputcsvdwell(YHZ68)
outputcsvdwell(YHZ75)

function datasummary = outputcsvdwell(dataset)
filename = sprintf('%s_fitting_results.csv','YHZ75');
datasummary = [];
for celliter = [1,2,3,4,5]
[beta1,CI_beta1] = cal_log_exp(dataset{celliter,1});
datasummary(celliter,1) = beta1;
datasummary(celliter,2) = CI_beta1(1);
datasummary(celliter,3) = CI_beta1(2);
end
csvwrite(filename,datasummary);
end


function out = getVarName(var)
    out = inputname(1);
end

function [beta1_e,CI_beta1_e] = cal_log_exp(datasets)

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
t_crit1 = tinv(0.95, df1);

CI_beta1 = [beta1 - t_crit1*SE_beta1, beta1 + t_crit1*SE_beta1];
CI_beta1_e = (CI_beta1);
beta1_e = (beta1);

end



