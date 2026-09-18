% cell size comparsion 
clear;close all;
celltypeidx = [3:1:7]; 
cellsize = cell(1,1);
for celltypeidxiter = 1:length(celltypeidx)
dirstring = sprintf('*merge_Darea_RG*-%d.mat',celltypeidx(celltypeidxiter));
filelist = dir(dirstring);
%%
resultsdata = cell(1,1);
for dirlength = 1:length(filelist)
    resultsfilename = filelist(dirlength).name;
    resultsdata{dirlength,1} = load(resultsfilename);
end

b = cell(1,1);
for i = 1:length(resultsdata)
    a = resultsdata{i,1};
    b{i,1} = a.Darea2;
end
c = vertcat(b{:});

cellsize{celltypeidxiter,1} = cell2mat(c(:,4));

end
%%
cell1 = cellsize{1,1};
cell2 = cellsize{2,1};
cell3 = cellsize{3,1};
cell4 = cellsize{4,1};
cell5 = cellsize{5,1};

x = [cellsize{1,1};cellsize{2,1};cellsize{3,1};cellsize{4,1};cellsize{5,1}];
g1 = repmat({'Mother-early'},length(cell1),1);
g2 = repmat({'Daughter-early'},length(cell2),1);
g3 = repmat({'Mother-late'},length(cell3),1);
g4 = repmat({'Daughter-late'},length(cell4),1);
g5 = repmat({'Individual'},length(cell5),1);

g = [g1; g2; g3; g4; g5];
boxplot(x(:,1),g);
%%
x1 = rand(5,1);
x2 = rand(10,1);
x3 = rand(15,1);
x = [x1; x2; x3];