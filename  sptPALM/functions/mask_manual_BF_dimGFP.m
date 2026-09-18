function [bf_label, bfPixelList] = mask_manual_BF_dimGFP(BF, BF_ori,savename)
close all;
blocations = bwboundaries(BF,'noholes');
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
Iblur3 = imbinarize(Iblur, 0.5);

bf_label = bwlabel(Iblur3,8);
stats = regionprops(bf_label,Iblur3,'Image','PixelList','PixelValues');
bfPixelValues = {stats.PixelValues}.';
bfPixelValues = {stats.PixelValues}.';
bfPixelList = {stats.PixelList}.';


save_cell = strrep(savename,'mask','bfPixelList');
save(save_cell, 'bfPixelList');
save_image = strrep(savename,'mask','manul_mask');
save(save_image, 'bf_label');
save_image = strrep(savename,'mask','freehand_mask');
save(save_image, 'hfhs');
%%
figure 
imshow(Iblur3);



