function [gfp_label,gfpPixelList] = mask_manual_GFP_dimGFP(GFPfilename,BF,savename,trackfile)
% read the GFP file 
close all;
tiff_stack = imread(GFPfilename, 1) ; % read in first image
%concatenate each successive tiff to tiff_stack
tiff_info = imfinfo(GFPfilename);

for ii = 2:length(tiff_info)
    temp_tiff = imread(GFPfilename, ii);
    tiff_stack = cat(3 , tiff_stack, temp_tiff);
end
mipGFP = max(tiff_stack,[],3);

figure
imagesc(mipGFP);
colormap
hold on

% read the BF file 
blocations = bwboundaries(BF,'noholes');

% read the track file 
data = readmatrix(trackfile{1});
data2 = sortrows(data,[1,2]);
spots_x = data2(:,3);
spots_y = data2(:,4);
% draw a circle centered at the BF cell;
for ind = 1:numel(blocations)
    xDim = 512;
    yDim = 512;
    pos = blocations{ind};
     % Convert to x,y order.
    pos = fliplr(pos);
    posx = pos(:,1);
    posy = pos(:,2);
centers(1) = mean(posx);
centers(2) = mean(posy);
radii = 1.35*(length(pos)/3.14).^0.5;
xc = centers(:,1);
yc = centers(:,2);
[xx,yy] = meshgrid(1:yDim,1:xDim);
gfpmask = false(xDim,yDim);
for ii = 1:numel(radii)
	gfpmask = gfpmask | hypot(xx - xc(ii), yy - yc(ii)) <= radii(ii);
end

   gfpblocations = bwboundaries(gfpmask,'noholes');
    gfppos = gfpblocations{1};
     % Convert to x,y order.
    gfppos = fliplr(gfppos);
    waypointsgfp = zeros(length(gfppos),1);
    waypointsgfp(round(length(waypointsgfp)*0.25)) = 0;
    waypointsgfp(round(length(waypointsgfp)*0.75)) = 0;
    drawfreehand('Position', gfppos,'Waypoints',logical(waypointsgfp),'Color','r');
end
hold on;
scatter(spots_x,spots_y,1,'r','MarkerEdgeAlpha',0.1); 
hold off;

%%
hfhs = findobj(gca, 'Type', 'images.roi.Freehand');
editedMask = false(size(BF));
for ind = 1:numel(hfhs)
    % Accumulate the mask from each ROI
    editedMask = editedMask | hfhs(ind).createMask();
    boundaryLocation = hfhs(ind).Position;
    boundaryLocationint = round(boundaryLocation);
    bInds = sub2ind(size(BF), boundaryLocationint(:,2), boundaryLocationint(:,1));
    sxeditedMask(bInds) = true;
end
editedMask2 = double(editedMask);
Iblur = imgaussfilt(editedMask2,1); 
Iblur3 = imbinarize(Iblur, 0.5);

gfp_label = bwlabel(Iblur3,8);
stats = regionprops(gfp_label,Iblur3,'Image','PixelList','PixelValues');
gfpPixelValues = {stats.PixelValues}.';
gfpPixelList = {stats.PixelList}.';


save_cell = strrep(savename,'BF','gfpPixelList');
save(save_cell, 'gfpPixelList');
save_image = strrep(savename,'BF','manul_mask_gfp');
save(save_image, 'gfp_label');
save_image = strrep(savename,'BF','freehand_mask_gfp');
save(save_image, 'hfhs');
%%
figure 
imshow(Iblur3);
end


