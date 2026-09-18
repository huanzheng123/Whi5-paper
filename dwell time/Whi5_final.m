% plot the track lenngth vs the cell volumn 
filename= '*merge_Darea_Dwellres*%d.mat';
savename = 'Whi5_bound_time_results.mat';

celltype = [3 5];
sizethre = [30 130];
[boundtime_extr1,cellvol_1] = boundtimecal(celltype,sizethre,filename);

celltype = [4 6];
sizethre = [10 60];
[boundtime_extr2,cellvol_2] = boundtimecal(celltype,sizethre,filename);

celltype = [4 6];
sizethre = [10 40];
boundtime_extr3 = boundtimecal(celltype,sizethre,filename);

celltype = [4 6];
sizethre = [40 60];
boundtime_extr4 = boundtimecal(celltype,sizethre,filename);

celltype = [7];
sizethre = [15 100];
[boundtime_extr5,cellvol_5] = boundtimecal(celltype,sizethre,filename);

figure
histogram(boundtime_extr3,'BinWidth',1,'Normalization','pdf','FaceColor','green');
hold on 
histogram(boundtime_extr4,'BinWidth',1,'Normalization','pdf','FaceColor','black');
%%
results = cell(1,1);
results{1,1} = boundtime_extr1;
results{2,1} = boundtime_extr2;
results{3,1} = boundtime_extr3;
results{4,1} = boundtime_extr4;
results{5,1} = boundtime_extr5;
results{1,2} = cellvol_1;
results{2,2} = cellvol_2;
results{5,2} = cellvol_5;

save (savename,'results');
%%
close all;

duration = [boundtime_extr1',boundtime_extr2',boundtime_extr3',boundtime_extr4'];
DURM = repmat({'Mother cell'},length(boundtime_extr1),1);
DURD = repmat({'Daughter cell'},length(boundtime_extr2),1);
DURdau = repmat({'Daughter cell (small)'},length(boundtime_extr3),1);
DURdau2 = repmat({'Daughter cell (large)'},length(boundtime_extr4),1);
DUR = vertcat(DURM,DURD,DURdau,DURdau2);
grouporder2={'Mother cell','Daughter cell','Daughter cell (small)','Daughter cell (large)'};
figure
violinplot(duration',DUR,'GroupOrder',grouporder2);


function [boundtime_extr,cellvol2] = boundtimecal(celltype,sizethre,filename)
binnum = length(sizethre) - 1;
% combine the results first;  
% prepare the data to analysis 
if length(celltype) == 1
    Dareafilename = sprintf(filename,celltype);
    Dareafile = dir(Dareafilename);
    Darea = cell(length(Dareafile),1);
for i = 1:length(Dareafile)
Darea{i} = importdata(Dareafile(i).name);
end
DareaAll = vertcat(Darea{:});

else
    for celltypeiter = 1:length(celltype)
        Dareafilename = sprintf(filename,celltype(celltypeiter));
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

% plot the bound time different 
boundtime_extr = extrinfo (dwellbinned{1,1});
boundtime_extr = boundtime_extr.*0.5; 
end

function boundtime_extr = extrinfo(strucutrelist)
dwellbinned2 = strucutrelist;
boundtimeall = {dwellbinned2(:).boundtime};
boundtime_extr = horzcat(boundtimeall{:})';

a = find(boundtime_extr >25); % delete the outliers 
boundtime_extr(a) = [];

end