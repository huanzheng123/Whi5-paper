% check BF seg 
clear; close all;
% switchs 
trackfilter = 1; % 1 means do the trackfilter 
BFGFPmanualmask = 1; % 1 means do the bf mask, otherwise load the cellarea file; 
DhistMinSteps = 5;
% input the data
fileall = [1]; % input the index of field of view 
cellpro = [7];
for fileiter = 1:length(fileall)
    for cellproiter = 1:length(cellpro)
        close all;
        fileidx = fileall(fileiter);
        celltypeiter = cellpro(cellproiter);
%function GFPseg()
if fileidx < 10
    BF_filename = sprintf('75_20220907_MS574_0%d_w2T2_merge_GFP-%d.mat',fileidx,celltypeiter);
    BF_file= imread(BF_filename);
    BF_file= BF_file(:,:,1);
    BF_oriname = sprintf('75_20220907_MS574_BF_0%d.tif',fileidx);
    BF_ori = imread(BF_oriname);

end
%%
% do the cell mask 
if BFGFPmanualmask == 1
[BFManuMask,bfPixelList] = mask_manual_BF_check(BF_file, BF_ori);
% do the nuclear mask
% mip of GFP, resize the hfh then manually mask 



end
    end
end

function [bf_label, bfPixelList] = mask_manual_BF_check(BF, BF_ori)
BF_file2 = imbinarize(BF);
BF_file3 = ~BF_file2;
BF_file4 = imfill(BF_file3,'holes');
blocations = bwboundaries(BF_file4,'noholes');
figure
imshow(BF_ori);
colormap(gray(256));
hold on
for ind = 1:numel(blocations)
    % Convert to x,y order.
    pos = blocations{ind};
    pos = fliplr(pos);
    % Create a freehand ROI.
    drawfreehand('Position', pos, 'Smoothing');
end
hold off;
% figure
% subplot(2,1,1);
% imagesc(BF_ori);
% colormap(gray(256));
% subplot(2,1,2);
% fuse = imfuse(GFP2,BF);
% a = imshow(fuse);
% waitfor(a);
%%
hfhs = findobj(gca, 'Type', 'images.roi.Freehand');
editedMask = false(size(BF));
for ind = 1:numel(hfhs)
    % Accumulate the mask from each ROI
    editedMask = editedMask | hfhs(ind).createMask();
    boundaryLocation = hfhs(ind).Position;
    boundaryLocationint = round(boundaryLocation);
    bInds = sub2ind(size(BF), boundaryLocationint(:,2), boundaryLocationint(:,1));
    editedMask(bInds) = true;
end
editedMask2 = double(editedMask);
Iblur = imgaussfilt(editedMask2,1); 
Iblur3 = imbinarize(editedMask2, 0.5);

bf_label = bwlabel(Iblur3,8);
stats = regionprops(bf_label,Iblur3,'Image','PixelList','PixelValues');
bfPixelValues = {stats.PixelValues}.';
bfPixelValues = {stats.PixelValues}.';
bfPixelList = {stats.PixelList}.';

%%
figure 
imagesc(Iblur3);
end



