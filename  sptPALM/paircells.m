clear; close all;
% pair the daughter/mother cells
% show 2 BF seg images and label in middle; overlap 2 images; 
field = [1:1:5];
pairchart = cell(4,length(field));
celltype = 1; % 1 means cell type 1/2, which are files with 3/4 suffix; 
% 3 means cell type 3/4, which are files with 5/6 suffix; 
for fielditer = 1:length(field)
    filename = '42_20201211_merge_manul_mask_Z_0%d-%d.mat';
cellgroup1 = sprintf(filename,field(fielditer),celltype+2);
cellgroup2 = sprintf(filename,field(fielditer),celltype+3);
pixelgroup1 = strrep(cellgroup1,'manul_mask','bfPixelList');
pixelgroup2 = strrep(cellgroup2,'manul_mask','bfPixelList');
bf_1 = importdata(cellgroup1);
bf_2 = importdata(cellgroup2);
pix_1 = importdata(pixelgroup1);
pix_2 = importdata(pixelgroup2);
c = imfuse(bf_1,bf_2);
imshow(c);
hold on
labelpos1 = [];
labelpos2 = [];
for i = 1:min(size(pix_1,1),size(pix_2,1))
labelpos1 = pix_1{i,1};
text(labelpos1(1,1),labelpos1(1,2),num2str(i),'Color','red','FontSize',14);
labelpos2 = pix_2{i,1};
text(labelpos2(1,1),labelpos2(1,2),num2str(i),'Color','yellow','FontSize',14);
hold on
end
a = 1;
pairchart{celltype,field(fielditer)} = [1:1:i];
pairchart{celltype+1,field(fielditer)} = [1:1:i];
end
%%
celltype = 3; % 1 means cell type 1/2, which are files with 3/4 suffix; 
% 3 means cell type 3/4, which are files with 5/6 suffix; 
for fielditer = 1:length(field)
cellgroup1 = sprintf(filename,field(fielditer),celltype+2);
cellgroup2 = sprintf(filename,field(fielditer),celltype+3);
pixelgroup1 = strrep(cellgroup1,'manul_mask','bfPixelList');
pixelgroup2 = strrep(cellgroup2,'manul_mask','bfPixelList');
bf_1 = importdata(cellgroup1);
bf_2 = importdata(cellgroup2);
pix_1 = importdata(pixelgroup1);
pix_2 = importdata(pixelgroup2);
c = imfuse(bf_1,bf_2);
imshow(c);
hold on
labelpos1 = [];
labelpos2 = [];
for i = 1:min(size(pix_1,1),size(pix_2,1))
labelpos1 = pix_1{i,1};
text(labelpos1(1,1),labelpos1(1,2),num2str(i),'Color','red','FontSize',14);
labelpos2 = pix_2{i,1};
text(labelpos2(1,1),labelpos2(1,2),num2str(i),'Color','yellow','FontSize',14);
hold on
end
a = 1;
pairchart{celltype,field(fielditer)} = [1:1:i];
pairchart{celltype+1,field(fielditer)} = [1:1:i];
end
