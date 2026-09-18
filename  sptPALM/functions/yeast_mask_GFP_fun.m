



function [GFP, nulPL] = yeast_mask_GFP_fun(filename)
imageread=imread(filename);
bw = imgaussfilt(imageread,1.5); %Could change it to 2
image_bi = imbinarize(bw,'adaptive');
imagecell = bwareaopen(image_bi,90);
image_label = bwlabel(imagecell,8);
stats = regionprops(image_label,imageread,'Image','PixelList','PixelValues');
PixelValues = {stats.PixelValues}.';
image_binary = cell(length(PixelValues),1);
PixelValues = {stats.PixelValues}.';
PixelList = {stats.PixelList}.';
Pixel_binary = cell(length(PixelList),1);
% save GFP_mask_01_filtered.mat image_label;
GFP = image_label;
nulPL = PixelList;
end
