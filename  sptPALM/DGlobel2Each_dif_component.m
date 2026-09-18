clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%% parameters and input the data  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dareafile = uigetfile('*Darea_RG*','MultiSelect','on');
% Trackfile = uigetfile('*Trackinfo*','MultiSelect','on');
% Dareafilename = sprintf('*merge_Darea_RG*.mat'); %,celltypeidx);  %global
% Trackfilename = sprintf('*merge_trackinfo_RG*.mat'); %,celltypeidx);
cellall = [3 4 5 6 7];
dataplot = cell(1,1);
for celliter = 1:length(cellall)
    celltypeidx = cellall(celliter);
Dareafilename = sprintf('*merge_Darea_RG*%d.mat',celltypeidx); % grouped
Trackfilename = sprintf('*merge_trackinfo_RG*%d.mat',celltypeidx);
Dareafile = dir(Dareafilename);
Trackfile = dir(Trackfilename);
%%
forcefit = 2;
GMMfix = 1 % decide if test different GMM model or not, if set to 0, use the default GMM model
% number to fit GMM
maxiteraterGMM = 4;
DefaultGMM = 3;
numcomp = 2;
gofthreshold = 0.5;
imagesave = 0;
boundhardthre = exp(0.1);
figuresdir = '/Users/zhenghuan/Downloads/P_20190404_FinalResult_mutipara/';
minnumconsi = 10; % threshold number of tracks in each cell of calculate the distribution of D
% dif_idx = 1; % input the group idx,1 means bound, 2 means slow diffusive, 3 means fast diffusive 
% Dareafile = uigetfile('*Darea*.mat', 'Multiselect', 'on');
% Trackfile = uigetfile('*TrackInfo*.mat', 'Multiselect', 'on');
% goffile = uigetfile('*gof*.mat', 'Multiselect', 'on');
savename = Dareafile.name;
%%%%%%%%%%%%%%%%%%%%%%%%%% calculation %%%%%%%%%%%%%%%%%%%%%%%%%%
% combine all data in all the field of view 
Darea = cell(length(Dareafile),1);
for i = 1:length(Dareafile)
Darea{i} = importdata(Dareafile(i).name);
end
DareaAll = vertcat(Darea{:});
%%
TrackInfo = cell(length(Trackfile),1);
for i = 1:length(Trackfile)
TrackInfo{i} = importdata(Trackfile(i).name);
end
TrackInfo2 = vertcat(TrackInfo{:});
TrackInfoSort = sortrows(TrackInfo2,4);

%%
% Calculate the bound proportion and index of group for tracks 
% DataInx: 1. Field of view number. 2. The index of BF cell correspond to the 
% saved labelled images from cell_filter.m script; 3. The index of GFP cell
% correspond to the index in track filter document (*tracks_record.mat); 4. 
% The area of cell; 5-7. the index in trackinfo file 
% 8. The D value of tracks inside the cell. 9. the index of classfication
% group
[DataInx] = GMM_log_set_ini_hard_thre(DareaAll,savename,DefaultGMM,boundhardthre);
DataInxvol{celliter,1} = unique(DataInx(:,1:5),'rows');
DataInx(:,5:9) = DataInx(:,6:10);
DataInx(:,10) = [];
%%
%find the bound component's index 
DMean = cell(numcomp,1);
for CluIndex = 1: numcomp % default value is 3 
    Clu = find (DataInx(:,9) == CluIndex);
    DVal = DataInx(:,8);
    CluVal = DVal(Clu);
    DMean{CluIndex} = mean(CluVal);
end
%%
DMean2 = cell2mat(DMean);
[DMeansort, DMeansortInx] = sort(DMean2);
%% 
% calculate the bound fraction in single cell
DataSort = sortrows(DataInx,[1 2]);
BoundData = cell(size(DataInx,1),3);
BoundDataIdx = 1;
FieVal = unique(DataInx(:,1));
for FieInx = FieVal(1):FieVal(length(FieVal))
    a = find(DataSort(:,1) == FieInx);
    DataField = DataSort(a,:);
    [CellVal,ia,ic] = unique(DataField(:,2));
    CellVal_counts = accumarray(ic,1);
    for CellNum = 1:length(CellVal_counts)
        if CellVal_counts(CellNum) <  minnumconsi
            continue
        else
            b = find(DataField(:,2) == CellVal(CellNum));
            DataCell = DataField(b,:);        
            for boundiditer = 1:numcomp
            BoundNum = find (DataCell(:,9) == DMeansortInx(boundiditer));
            BoundFra(boundiditer) = length(BoundNum)/CellVal_counts(CellNum);
            end
            BoundData{BoundDataIdx,1} = FieInx;
            BoundData{BoundDataIdx,2} = CellVal(CellNum);
            BoundData{BoundDataIdx,3} = BoundFra;
            BoundDataIdx = BoundDataIdx+1;
        end
    end
end
BoundData(all(cellfun('isempty',BoundData),2),:) = [];
%%
% compact the data, add the bound fraction in to the result
DataSort2 = cell(size(DataSort,1),size(DataSort,2)-2);
DS2i = 1;
for ComFC = 1: size(BoundData,1)
    CellArea = find(DataSort(:,1)==BoundData{ComFC,1} & DataSort(:,2)==BoundData{ComFC,2});
    DataSort2(DS2i,1:4) = num2cell(DataSort(CellArea(1),1:4));
    DDataSort = DataSort(CellArea,8);
    TrackDataSort = DataSort(CellArea,5:7);
    DataSort2{DS2i,5} = [DDataSort TrackDataSort];
    DataSort2{DS2i,6} = DataSort(CellArea,9);
    DataSort2{DS2i,7} = BoundData{ComFC,3};
    DS2i = DS2i + 1;
end
DataSort2 = DataSort2(~any(cellfun('isempty',DataSort2), 2), :);

%%
% save the final result 
Result = struct();
settings = struct();
settings.DefaultGMM = DefaultGMM; 
settings.minnumconsi = minnumconsi;
data = struct();
data.InputData = Dareafile; 
data.Result = DataSort2;
% DataSort2: 1. Field of view number. 2. The index of BF cell correspond to the 
% saved labelled images from cell_filter.m script; 3. The index of GFP cell
% correspond to the index in track filter document (*tracks_record.mat); 4. 
% The area of cell; 5. Matrix, the D value, the index in trackinfo file 
% 6. the index of classfication. 7. the bound proportion in each cell 
% data.BoundIndex = BoundIdx;
% BoundIdx: the bound index in GMM model 
data.TrackInfo = TrackInfoSort;
% TrackInfoSort: the trackinfo, the field of view number 
Result.settings = settings;
Result.data = data;

savenameWS = strrep(savename,'Darea','WorkSpace_GroupTrend');
savenameFR = strrep(savename,'Darea','FinalResult_GroupTrend');
mkdir('Final_Result');
path = [pwd, filesep, 'Final_Result'];
matfile = fullfile(path, savenameFR);
save(matfile,'Result');
matfile2 = fullfile(path, savenameWS);
save(matfile2);
%%
%plot the final results 
xcellarea = cell2mat(DataSort2(:,4))/100;
yboundfra = cell2mat(DataSort2(:,7));
y2 = yboundfra(:,1)+ yboundfra(:,2);
ysum = yboundfra(:,1);

dataplot{celliter,1} = xcellarea; 
dataplot{celliter,2} = ysum;
dataplot{celliter,3} = DataSort2(:,1:4);
dataplot{celliter,4} = cellfun('length',DataSort2(:,6));
dataplot{celliter,5} = DataInx;
dataplot(:,6) = DataInxvol;

for figureid = 1:size(ysum,2)
figure
plot(xcellarea,ysum(:,figureid),'r.','MarkerSize',17);
title('The-G proportion of bound verse cell area');
xlabel('cell area(um^2)');
ylabel('Bound proportion');
axis([min(cell2mat(DataSort2(:,4))/100) max(cell2mat(DataSort2(:,4))/100) 0 1]);
set(gca,'FontSize',17);
if imagesave
    strname = sprintf('Figure_proportion_#%d',figureid);
    savenameFG = strrep(savename,'Darea',strname);
    savenameFG_2 = strrep(savenameFG,'.mat','.jpg');
    saveas(gcf,strcat(figuresdir, savenameFG_2), 'jpeg');
end
end
end
savenameplot = strrep(savename,'Darea','Data2Plot');
save(savenameplot,'dataplot');