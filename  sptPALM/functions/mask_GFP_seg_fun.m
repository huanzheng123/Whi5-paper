% GFP channel nuclear mask 
function [gfp_label,gfpPixelList] = mask_GFP_seg_fun(GFPfilename,savename)
close all;
tiff_stack = imread(GFPfilename, 1) ; % read in first image
%concatenate each successive tiff to tiff_stack
tiff_info = imfinfo(GFPfilename);

for ii = 2:length(tiff_info)
    temp_tiff = imread(GFPfilename, ii);
    tiff_stack = cat(3 , tiff_stack, temp_tiff);
end
mipGFP = max(tiff_stack,[],3);
bw = imgaussfilt(mipGFP,1.5); %Could change it to 2
image_bi = imbinarize(bw,'adaptive',Sensitivity=0.3);
image_label = bwlabel(image_bi,8);

stats = regionprops(image_label,mipGFP,'Image','PixelList','PixelValues');
PixelValues = {stats.PixelValues}.';
image_binary = cell(length(PixelValues),1);
PixelValues = {stats.PixelValues}.';
PixelList = {stats.PixelList}.';
Pixel_binary = cell(length(PixelList),1);
% save GFP_mask_01_filtered.mat image_label;
gfp_label = image_label;
gfpPixelList = PixelList;
% save(savename,"GFP");
imshow(gfp_label);
pause(1);
save_cell = strrep(savename,'BF','gfpPixelList');
save(save_cell, 'gfpPixelList');
save_image = strrep(savename,'BF','manul_mask_gfp');
save(save_image, 'gfp_label');
end