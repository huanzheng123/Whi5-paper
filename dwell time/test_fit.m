% to do the stat test 
% test the linear fit and compare the significancy

[f,x] = ecdf((boundtime_extr1),'Function','survivor','Alpha',0.01,'Bounds','on');
aaa = find(f == 0);
x(aaa) = [];
f(aaa) = [];

f2 = log(f);
aaa = find(f2 == 0);
x(aaa) = [];
f2(aaa) = [];

aaa = find(x>15);
x(aaa) = [];
f2(aaa) = [];

mdl1 = fitlm(x, f2);


%%
[fa,xa] = ecdf((boundtime_extr3),'Function','survivor','Alpha',0.01,'Bounds','on');
aaa = find(fa == 0);
xa(aaa) = [];
fa(aaa) = [];

f2a = log(fa);
aaa = find(f2a == 0);
xa(aaa) = [];
f2a(aaa) = [];

aaa = find(xa>15);
xa(aaa) = [];
f2a(aaa) = [];

mdl2 = fitlm(xa, f2a);

beta1 = mdl1.Coefficients.Estimate(2);
SE_beta1 = mdl1.Coefficients.SE(2);

beta2 = mdl2.Coefficients.Estimate(2);
SE_beta2 = mdl2.Coefficients.SE(2);

df1 = length(f2) - 2;
df2 = length(f2a) - 2;

t_crit1 = tinv(0.95, df1);
t_crit2 = tinv(0.95, df2);

CI_beta1 = [beta1 - t_crit1*SE_beta1, beta1 + t_crit1*SE_beta1];
CI_beta2 = [beta2 - t_crit2*SE_beta2, beta2 + t_crit2*SE_beta2];

fprintf('95%% CI for slope of first regression: [%f, %f]\n', CI_beta1(1), CI_beta1(2));
fprintf('95%% CI for slope of second regression: [%f, %f]\n', CI_beta2(1), CI_beta2(2));
