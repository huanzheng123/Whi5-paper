function [DataInx,GMMboundfraction,GMMboundmu] = GMM_log_set_ini_hard_thre(data,filename,numcomp,hardthre)
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
Mu = [-0.33;0.67];
Sigma(:,:,1) = 0.25;
Sigma(:,:,2) = 0.5;
PComponents = [0.5,0.5];
S = struct('mu',Mu,'Sigma',Sigma,'ComponentProportion',PComponents);
%
options = statset('Display','off','MaxIter',1000,'UseParallel',true);
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
GMMboundfraction = BestModel.ComponentProportion;
GMMboundmu = BestModel.mu;
%%
% % % plot the fitting
DMean = cell(numcomp,1);
for CluIndex = 1: numcomp % default value is 3 
    Clu = find (DataInx(:,10) == CluIndex);
    DVal = DataInx(:,9);
    CluVal = DVal(Clu);
    DMean{CluIndex} = mean(CluVal);
end

DMean2 = cell2mat(DMean);
[DMeansort, DMeansortInx] = sort(DMean2);

% plot the histogram 
BoundNumAll = find (DataInx(:,10) == DMeansortInx(1));
boundplot = DataInx(BoundNumAll,:);

BoundNumAll2 = find (DataInx(:,10) == DMeansortInx(2));
boundplot2 = DataInx(BoundNumAll2,:);

dataplot = log(boundplot(:,9));
dataplot2 = log(boundplot2(:,9));
figure
% histogram (dataplot,'BinWidth',0.1,'BinLimits',[-3,3],'FaceColor','blue');
% hold on
% histogram (dataplot2,'BinWidth',0.1,'BinLimits',[-3,3],'FaceColor','blue');
fitplotdata = log(DataInx(:,9));

histogram (fitplotdata,'BinWidth',0.1,'Normalization','pdf','FaceColor','blue');

boundfratext = BestModel.ComponentProportion;
boundfratext = boundfratext(1);
mutext = BestModel.mu;
mu1 = mutext(1);
mu2 = mutext(2);

sigmatext = BestModel.Sigma;
sigmatext1 = sigmatext(1);
sigmatext2 = sigmatext(2);

prportiontext = BestModel.ComponentProportion;
prportiontext1 = prportiontext(1);
prportiontext2 = prportiontext(2);
addtext = sprintf('mu % .4f,% .4f; fraction % .4f', mu1,mu2,boundfratext);

ylabel('pdf','FontSize',17);
xlabel(addtext,'FontSize',17);
% 
hold on
x= [min(fitplotdata):0.01:max(fitplotdata)];
x = x(:);
y = pdf(BestModel, x);
plot(x,y,'LineWidth',2);
hold off;

hold on 
gm = gmdistribution(mu1,sigmatext1);
y = prportiontext1*pdf(gm,x);
plot(x,y,'LineWidth',2);
hold off;

hold on 
gm = gmdistribution(mu2,sigmatext2);
y = prportiontext2*pdf(gm,x);
plot(x,y,'LineWidth',2);
hold off;

%%
% save the bestmodel and fitting plots and the data used in fitting
savenameModel = strrep(filename,'Darea','GMMAll_binned');
savenameData = strrep(filename,'Darea','DataAll_binned');
mkdir('GMM_Model');
path = [pwd, filesep, 'GMM_Model'];
matfile = fullfile(path, savenameModel);
matfile2 = fullfile(path, savenameData);
save(matfile,'BestModel');
save(matfile2, 'DataInx');
%%
    strname = sprintf('GMM_fit');
    savenameFG = strrep(filename,'Darea',strname);
    savenameFG_2 = strrep(savenameFG,'.mat','.jpg');
    matfile3 = fullfile(path, savenameFG_2);
    saveas(gcf,matfile3, 'jpeg');

    %%
    for j = 1:size(DataInx,1)
    if DataInx(j,9) < hardthre
        DataInx(j,10) = 1;
    else
        DataInx(j,10) = 2;
    end
    end
