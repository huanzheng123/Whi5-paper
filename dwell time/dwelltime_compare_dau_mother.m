%dwelltime_compare_dau_mother
clear; close all;
celltype = [2,4];
sizethre = [10,50,50];
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

boundtime_extr = extrinfo (dwellbinned{1,1});
%%
celltype = [1,3];
sizethre = [30,80,80];
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

boundtime_extr2 = extrinfo (dwellbinned{1,1});
%%
boundtime_extr = boundtime_extr.*0.5;
boundtime_extr2 = boundtime_extr2.*0.5;
%%
duration = [boundtime_extr',boundtime_extr2'];
DURM = repmat({'Daughter cell (10fL -- 30fL)'},length(boundtime_extr),1);
DURD = repmat({'Daughter cell (30fL -- 50fL)'},length(boundtime_extr2),1);
DUR = vertcat(DURM,DURD);
grouporder2={'Daughter cell (10fL -- 30fL)','Daughter cell (30fL -- 50fL)'};
figure
violinplot(duration',DUR,'GroupOrder',grouporder2);

function boundtime_extr = extrinfo(strucutrelist)
dwellbinned2 = strucutrelist;
boundtimeall = {dwellbinned2(:).boundtime};
boundtime_extr = horzcat(boundtimeall{:})';
end