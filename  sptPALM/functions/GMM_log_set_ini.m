function [DataInx,numcomp] = GMM_log_set_ini (data,forcefit,filename,figuresdir,imagesave)
% set localparameter
%GMM_models_tested = 4;
% unwind the cell array 
% DataInx: 1. Field of view number. 2. The index of BF cell correspond to the 
% saved labelled images from cell_filter.m script; 3. The index of GFP cell
% correspond to the index in track filter document (*tracks_record.mat); 4. 
% The area of cell; 5-7. the index in trackinfo file 
% 8. The D value of tracks inside the cell. 9. the index of classfication
% group
datarep = data(:,1:4);
datatrack = data(:,5);
cellrep = cell(size(datarep,1),1);
for i = 1:size(datarep,1)
    cellrep{i} = repmat(cell2mat(datarep(i,:)),size(datatrack{i},1),1);
end
%%
cellrep2 = vertcat(cellrep{:});
datatrack2 = vertcat(datatrack{:});
datatrack3 = datatrack2(:,1);
trackinfo = datatrack2(:,2:4);
cellrep3 = [cellrep2 trackinfo];
%%
% input the D
data_neg = find(datatrack3<=0); % del negetive data
datatrack3(data_neg) = [];
data_log = log(datatrack3);
% use BIC find the besxt fit GMM model
Mu = [-0.18;0.75];
Sigma(:,:,1) = 0.25;
Sigma(:,:,2) = 0.5;
PComponents = [0.5,0.5];
S = struct('mu',Mu,'Sigma',Sigma,'ComponentProportion',PComponents);
%
options = statset('Display','off','MaxIter',1000000,'UseParallel',true);
%GMModels = cell(1,GMM_models_tested);
%BIC = zeros(1,GMM_models_tested);
%%
%
%for k = 1:GMM_models_tested
    GMModels = fitgmdist(data_log,2,'Start',S,'Options',options,'Replicates',1,'SharedCovariance', false);
    BIC= GMModels.BIC;
%end
%[minBIC,numComponents_BIC] = min(BIC);
BestModel = GMModels;

%% 
% Return the bound d and fraction
[bound_log, bound_index] = min(BestModel.mu);
bound_d = exp(bound_log);
boundfra = BestModel.ComponentProportion(bound_index);
% return the index of the tracks 
[idx,nlogL] = cluster(BestModel,data_log);
% create the dataset 
DataInx = [cellrep3, datatrack3,idx];
numcomp = BestModel.NumComponents;
%%
% % % plot the fitting
data = data_log;
figure
histogram (data,'BinWidth',0.1,'Normalization','pdf','BinLimits',[-3,3]);
% hold on
% line([-2.2584 -2.2584],[0 0.3],'LineWidth',2,'Color','r');
%title('','FontSize',17);
ylabel('pdf','FontSize',17);
xlabel('D (a.u.)','FontSize',17);
% hold on
% x= [min(data):0.01:max(data)];
% x = x(:);
% y = pdf(BestModel, x);
% plot(x,y,'LineWidth',2);
% hold off;
%%
% save the bestmodel and fitting plots and the data used in fitting
savenameModel = strrep(filename,'Darea','GMMAll_slow_dif');
savenameData = strrep(filename,'Darea','DataAll_slow_dif');
mkdir('GMM_Model');
path = [pwd, filesep, 'GMM_Model'];
matfile = fullfile(path, savenameModel);
matfile2 = fullfile(path, savenameData);
save(matfile,'BestModel');
save(matfile2, 'DataInx');
%%
if imagesave
    strname = sprintf('GMM_fit');
    savenameFG = strrep(filename,'Darea',strname);
    savenameFG_2 = strrep(savenameFG,'.mat','.jpg');
    saveas(gcf,strcat(figuresdir, savenameFG_2), 'jpeg');
end
