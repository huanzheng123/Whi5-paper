% new verison of script to analysis whi5 nuclear marker spt dat 
clear; close all;
% switchs 
trackfilter = 1; % 1 means do the trackfilter 
BFGFPmanualmask = 1; % 1 means do the bf mask, otherwise load the cellarea file; 
DhistMinSteps = 4;
% input the data
fileall = [1 2 3 4]; % input the index of field of view 
cellpro = [4]; %
for fileiter = 1:length(fileall)
    for cellproiter = 1:length(cellpro)
        close all;
        fileidx = fileall(fileiter);
        celltypeiter = cellpro(cellproiter);
%function GFPseg()
if fileidx < 10
    BF_filename = sprintf('dwell_33_fake_20211015_MS162_merge_GFP_0%d_0%d.mat',fileidx,celltypeiter);
    BF_file= imread(BF_filename);
    BF_file= BF_file(:,:,1);
    BF_oriname = sprintf('dwell_33_20211021_MS162_0%d_w3BF.TIF',fileidx);
    BF_ori = imread(BF_oriname);
    savename = strrep(BF_filename,'GFP','BF');
    GFPfilename = sprintf('dwell_33_20211021_MS162_GFP_0%d.TIF',fileidx);
    trackfile{1} = sprintf('dwell_33_20211021_MS162_stream_0%d.TIF_spots.csv',fileidx);
%     trackfile{2} = strrep(trackfile{1},'.TIF','.tif');
%     trackfile{3} = strrep(trackfile{1},'.TIF','.tif');
%     trackfile{4} = strrep(trackfile{1},'.TIF','.tif');
else
    BF_filename = sprintf('20200615_ms564_bf_mask_Z Series_%d.mat',fileidx);
    BF_file = load(BF_filename);
    BF = BF_file.final_image;
    BF_oriname = sprintf('20200615_ms564_bf_%d.tif',fileidx);
    BF_ori = imread(BF_oriname);
    savename = BF_filename;
    GFPfilename = sprintf('20200615_ms564_gfp_Z Series_%d.tif',fileidx);
    trackfile{1} = sprintf('20200619_MS563_STREAM_thre220_0%d.tif_spots.csv',fileidx);
    trackfile{2} = sprintf('20200619_MS563_STREAM_thre220_0%d-file002.tif_spots.csv',fileidx);
    trackfile{3} = sprintf('20200619_MS563_STREAM_thre220_0%d-file003.tif_spots.csv',fileidx);
end
%%
% do the cell mask 
if BFGFPmanualmask == 1
[BFManuMask,bfPixelList] = mask_manual_BF_newSeg(BF_file,BF_ori,savename);
% do the nuclear mask
% mip of GFP, resize the hfh then manually mask 
%%
[GFPManuMask,gfpPixelList] = mask_manual_GFP_dimGFP(GFPfilename,BFManuMask,savename,trackfile);
%%
%%
% make cellarea file 
cellarea = zeros(length(gfpPixelList),3);
 for i = 1:length(gfpPixelList)
     nulpos = gfpPixelList{i,1};
     med_x = int16(median(nulpos(:,1)));
     nul_yid = find(nulpos(:,1)==med_x);
     nul_y = nulpos(:,2);
     med_y = int16(median(nul_y(nul_yid)));
     bf_index = BFManuMask(med_y,med_x);
     cellarea(i,2) = GFPManuMask(med_y,med_x);
     if bf_index~=0 
         cellarea(i,3) = length(bfPixelList{bf_index,1});
         cellarea(i,1) = bf_index;
     end
 end
 cellidxfind = find(cellarea(:,1) == 0);
 cellarea(cellidxfind,2) = 0;
 save_cellarea = strrep(savename,'BF','cell_area');
 save(save_cellarea, 'cellarea');
else
    save_cellarea = strrep(savename,'BF','cell_area');
    cellarea = load(save_cellarea);
end
    end
end
