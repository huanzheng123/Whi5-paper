% from the cell size and dwell time results, get binned dwell time results;
clear; close all;
celltype = [2 4];
sizethre = [10,35,60];
binnum = length(sizethre) - 1;
% combine the results first;  
% prepare the data to analysis 
if length(celltype) == 1
    Dareafilename = sprintf('*merge_Darea_Dwellres*%d.mat',celltype);
    Dareafile = dir(Dareafilename);
    Darea = cell(length(Dareafile),1);
for i = 1:length(Dareafile)
Darea{i} = importdata(Dareafile(i).name);
end
DareaAll = vertcat(Darea{:});

else
    for celltypeiter = 1:length(celltype)
        Dareafilename = sprintf('*merge_Darea_Dwellres*%d.mat',celltype(celltypeiter));
        Dareafile{celltypeiter,1} = dir(Dareafilename);
    end
    Dareafile2 = vertcat(Dareafile{:});
    Darea = cell(length(Dareafile2),1);
for i = 1:length(Dareafile2)
Darea{i} = importdata(Dareafile2(i).name);
end
DareaAll = vertcat(Darea{:});
end

DareaAllsize = DareaAll(:,4);
B = cellfun(@(x) x(:,2),DareaAllsize,'uni',0); % delete the cell area, only reserve cell volumn
cellvol = cellfun(@(x) x.*0.1^3,B,'uni',0); % convert to the proper unit fL;  
cellvol2 = cell2mat(cellvol); 

DareaAll2 = DareaAll; 
DareaAll2(:,4) = cellvol;

figure 
histogram(cellvol2,'BinWidth',3,'Normalization','pdf','FaceColor','blue');

% from the cell size threshold to extract the dwell time info 
DwellAllbinned = cell(binnum,1);
dwellbinned = cell(binnum,1);
for binnedgroupidx = 1:binnum
    DareaAllbinedidx = find(cellvol2 <= sizethre(binnedgroupidx+1) & cellvol2 > sizethre(binnedgroupidx));
    Dareabin = DareaAll2(DareaAllbinedidx,:);
    DareaAllbinned{binnedgroupidx,1} = Dareabin; 
    dwellbinned{binnedgroupidx,1} = vertcat(Dareabin{:,6});
end

% plot the bound time different 
boundtime_extr = extrinfo (dwellbinned{1,1});
boundtime_extr2 = extrinfo (dwellbinned{2,1});
figure
histogram(boundtime_extr,'BinWidth',1,'Normalization','pdf','FaceColor','green');
hold on 
histogram(boundtime_extr2,'BinWidth',1,'Normalization','pdf','FaceColor','black');
%%
boundtime_extr = boundtime_extr.*0.5;
boundtime_extr2 = boundtime_extr2.*0.5;
%%
celltype = [3 4 5 6 7];
sizethre = [];
binnum = length(sizethre) - 1;
% combine the results first;  
% prepare the data to analysis 
if length(celltype) == 1
    Dareafilename = sprintf('*merge_Darea_Dwellres*%d.mat',celltype);
    Dareafile = dir(Dareafilename);
    Darea = cell(length(Dareafile),1);
for i = 1:length(Dareafile)
Darea{i} = importdata(Dareafile(i).name);
end
DareaAll = vertcat(Darea{:});

else
    for celltypeiter = 1:length(celltype)
        Dareafilename = sprintf('*merge_Darea_Dwellres*%d.mat',celltype(celltypeiter));
        Dareafile{celltypeiter,1} = dir(Dareafilename);
    end
    Dareafile2 = vertcat(Dareafile{:});
    Darea = cell(length(Dareafile2),1);
for i = 1:length(Dareafile2)
Darea{i} = importdata(Dareafile2(i).name);
end
DareaAll = vertcat(Darea{:});
end

DareaAllsize = DareaAll(:,4);
B = cellfun(@(x) x(:,2),DareaAllsize,'uni',0); % delete the cell area, only reserve cell volumn
cellvol = cellfun(@(x) x.*0.10^3,B,'uni',0); % convert to the proper unit fL;  
cellvol2 = cell2mat(cellvol); 

DareaAll2 = DareaAll; 
DareaAll2(:,4) = cellvol;

figure 
histogram(cellvol2,'BinWidth',3,'Normalization','pdf','FaceColor','blue');
%%
% from the cell size threshold to extract the dwell time info 
DwellAllbinned = cell(binnum,1);
dwellbinned = cell(binnum,1);
for binnedgroupidx = 1:binnum
    DareaAllbinedidx = find(cellvol2 <= sizethre(binnedgroupidx+1) & cellvol2 > sizethre(binnedgroupidx));
    Dareabin = DareaAll2(DareaAllbinedidx,:);
    DareaAllbinned{binnedgroupidx,1} = Dareabin; 
    dwellbinned{binnedgroupidx,1} = vertcat(Dareabin{:,6});
end

boundtime_extr3 = extrinfo (dwellbinned{1,1});
boundtime_extr3 = boundtime_extr3.*0.5;


%%
%figure 2: track duration survival plot
duration = [boundtime_extr',boundtime_extr2',boundtime_extr3'];
DURM = repmat({'Daughter cell (10fL -- 30fL)'},length(boundtime_extr),1);
DURD = repmat({'Daughter cell (30fL -- 50fL)'},length(boundtime_extr2),1);
DURMother = repmat({'Mother cell'},length(boundtime_extr3),1);
DUR = vertcat(DURM,DURD,DURMother);
grouporder2={'Daughter cell (10fL -- 30fL)','Daughter cell (30fL -- 50fL)','Mother cell'};
figure
violinplot(duration',DUR,'GroupOrder',grouporder2);

%%
function boundtime_extr = extrinfo(strucutrelist)
dwellbinned2 = strucutrelist;
boundtimeall = {dwellbinned2(:).boundtime};
boundtime_extr = horzcat(boundtimeall{:})';
end













