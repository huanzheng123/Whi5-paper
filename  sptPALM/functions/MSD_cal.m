% MSD method
function [track_D] = MSD_cal(data, MethInd)
%input trajectories 
if nargin<2
  MethInd = 3; % least square fit of MSD 
end

if MethInd == 1 % use half of the tau value to calculate MSD 
%%
max_tau = round(length(data)/2);
         msd = zeros(max_tau-1,1);
         for i = 1:max_tau-1
             msd_sing = 0;
             for j = 1:length(data)-i
                 msd_sing = msd_sing + dist(data(j,1),data(j,2),data(j+i,1),data(j+i,2)).^2;
             end
             msd(i) = msd_sing/(length(data)-i);
         end
         a = max_tau-1;
         tau = [1:a]';
         
%         tau_inv = max_tau-tau;
% wei_sum = sum(tau_inv);
% Weights = tau_inv./wei_sum;
         
         MSD_fit = fittype('4*D*t^a + 4*b^2', 'dependent',{'y'},'independent',{'t'},'coefficients',{'D','a','b'});
         start = [D1, 1 loci];
         lower = [0, 0 locl];
         upper = [50, 2 locu];
         [MSD_params,gof,output] = fit(tau, msd, MSD_fit,'Start',start,'Lower',lower,'Upper',upper);
         fit_params = [MSD_params.D, MSD_params.a, MSD_params.b];
         tracks_MSDs{1,1} = msd;
         tracks_MSDs{1,2} = data;
         tracks_MSDs{1,3} = fit_params;
         tracks_MSDs{1,4} = gof;
         track_D = MSD_params.D;
% apply the radius gyration to make the classification 
elseif MethInd == 2
    data_mean = mean(data);
    sum_rg = 0;
    for i = 1:length(data)
        sum_rg = sum_rg + (data(i,1)-data_mean(1,1))^2 + (data(i,2)-data_mean(1,2))^2;
    end
    track_D = (sum_rg/i)^0.5;

% apply the weighted MSD to make the classification 
 elseif MethInd == 3
     max_tau = length(data);
msd = zeros(max_tau-1,1);
%%
for i = 1:max_tau-1
    msd_sing = 0;
    for j = 1:length(data)-i
        msd_sing = msd_sing + dist(data(j,1),data(j,2),data(j+i,1),data(j+i,2)).^2;
    end
    msd(i) = msd_sing/(length(data)-i);
end
%%
a = max_tau-1;
tau = [1:a]';

tau_inv = max_tau-tau;
wei_sum = sum(tau_inv);
Weights = tau_inv./wei_sum;
%%
   MSD_fit = fittype('4*D*t^a + 4*b^2', 'dependent',{'y'},'independent',{'t'},'coefficients',{'D','a','b'});
 
        start = [D1, 1 loci];
         lower = [0, 0 locl];
         upper = [50, 2 locu];
        [MSD_params,gof,output] = fit(tau, msd, MSD_fit,'Start',start,'Lower',lower,'Upper',upper,'Weights',Weights);
        fit_params = [MSD_params.D, MSD_params.a, MSD_params.b];
        %%
        %tracks_MSDs = struct('MSD_values', MSD, 'GoF',gof,'Parameter_Values', fit_params); 
        tracks_MSDs{1,1} = msd;
         tracks_MSDs{1,2} = data;
         tracks_MSDs{1,3} = fit_params;
         tracks_MSDs{1,4} = gof;
         track_D = MSD_params.D;

end
end


