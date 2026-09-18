% modify the copy number script to extract the 
% cell intensity extraction -- using the cylincal method 
clear, close all;
ra = 0.8798; % use this ra for now, need to estimate ra from MS10 later; 
rg = 0.1065; %from pure mNeonGreen data 
ratiov = 0.0995; % ratio of cell volume and nuclear volume in BY4741; 
pixelsize = 0.0725; % um per pixel; 

% input data: BF mask and GFP z a
% load the different BF mask layer
fileall = [1:1:50];
cellpro = [3 4 7]; 

for fileiter = 1:length(fileall)
    for cellproiter = 1:length(cellpro)
        close all;
        fileidx = fileall(fileiter);
        celltypeiter = cellpro(cellproiter);

if fileidx < 10
    BF_filename = sprintf('20230630_copynum_swi4_mNG_0%d_w2merge_GFP_cp_masks_%d.mat',fileidx,celltypeiter); % merged BF files
    
    GFP_filename_seg  = sprintf('20230630_copynum_swi4_mNG_0%d_w1_Cropped_MAX_GFP-488-561-610-75_cp_masks.tif',fileidx); % GFP segmentation image 
    GFP_seg = imread(GFP_filename_seg);
    savename = strrep(BF_filename,'GFP','BF');
    GFPfilename = sprintf('20230630_copynum_swi4_mNG_0%d_w1_Cropped_GFP-488-561-610-75new.tif',fileidx); % GFP cropped image
    RFPfilename = strrep(GFPfilename,'GFP','RFP'); 

    BFManuMask = sprintf("20230630_copynum_swi4_mNG_0%d_w2merge_manul_mask_cp_masks_%d.mat",fileidx,celltypeiter);
    BFManuMask = load(BFManuMask).bf_label;
    bfPixelList = sprintf("20230630_copynum_swi4_mNG_0%d_w2merge_bfPixelList_cp_masks_%d.mat",fileidx,celltypeiter);
    bfPixelList = load(bfPixelList).bfPixelList;

else
    BF_filename = sprintf('20230630_copynum_swi4_mNG_%d_w2merge_GFP_cp_masks_%d.mat',fileidx,celltypeiter);
   
    GFP_filename_seg  = sprintf('20230630_copynum_swi4_mNG_%d_w1_Cropped_MAX_GFP-488-561-610-75_cp_masks.tif',fileidx);
    GFP_seg = imread(GFP_filename_seg);
    savename = strrep(BF_filename,'GFP','BF');
    GFPfilename = sprintf('20230630_copynum_swi4_mNG_%d_w1_Cropped_GFP-488-561-610-75new.tif',fileidx);
    RFPfilename = strrep(GFPfilename,'GFP','RFP');

    BFManuMask = sprintf("20230630_copynum_swi4_mNG_%d_w2merge_manul_mask_cp_masks_%d.mat",fileidx,celltypeiter);
    BFManuMask = load(BFManuMask).bf_label;
    bfPixelList = sprintf("20230630_copynum_swi4_mNG_%d_w2merge_bfPixelList_cp_masks_%d.mat",fileidx,celltypeiter);
    bfPixelList = load(bfPixelList).bfPixelList;

end

    % get the BF and GFP seg

 
  [GFPManuMask,gfpPixelList,tiff_stack] = mask_GFP_copynum2(GFPfilename,GFP_seg,savename);
if max(max(BFManuMask)) == 0
    volume = [];
else
  volume = volume_estimate(BFManuMask);
end
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
         cellsizepixel = length(bfPixelList{bf_index,1});
         cellarea(i,3) = cellsizepixel.*pixelsize^2; 
         BFvolume = volume(bf_index);
         cellarea(i,4) = BFvolume.* pixelsize^3; 
         cellarea(i,1) = bf_index;
         nucsizepixel = length(gfpPixelList{i,1});
         cellarea(i,5) = nucsizepixel.*pixelsize^2; 
     end
 end

% Delete the zero cell area
Darea2 = cellarea;
areaind = find(Darea2(:,3)==0 | Darea2(:,2)==0);
Darea2(areaind,:)=[];

%%
% get the slice number and intensity for inquired GFP cell
for GFPiter = 1:size(Darea2,1)
    bfnum = Darea2(GFPiter,2);
    norintenGFPall = zlayerinten2(tiff_stack, RFPfilename, gfpPixelList,bfnum,ra,rg); 
    Darea2(GFPiter,6) = 0;
    Darea2(GFPiter,7) = 0;
    Darea2(GFPiter,8) = 0; 
    Darea2(GFPiter,9) = norintenGFPall;
    Darea2(GFPiter,10) = (norintenGFPall+201.96)./1742.1; % y = 1742.1X-201.96 mNG copy num relation with intensity;
    Darea2(GFPiter,11) = Darea2(GFPiter,10)./Darea2(GFPiter,4);
end

% save cell area: BF labeled mask and GFP labeled mask; 1.BF idx 2.nuclear idx 3. BF
% cell area 4. BF cell v 5. nuclear area 6.max intensity slice num 7. intensity of
% that slice 8.normailized intensity 9.estimated whole v intensity; 10.
% normalized copy number 11. concentration 

savename2 = strrep(savename,'BF','intensity_GFP_nuclear_cell_cylincal');

save(savename2,'Darea2');


    end
end

% structure of the data, BF idx, BF mask, BF cell size and volume, GFP cell
% idx, nuclear size and intensity, normalized intensity 

%%
% input coordinates and image, output intensity 
function norintenGFPall2 = zlayerinten2(tiff_stack, RFPfilename, bfPixelList,bfnum,ra,rg)
% read RFP channel stack 
tiff_stackRFP = imread(RFPfilename, 1) ; % read in first image
%concatenate each successive tiff to tiff_stack
tiff_infoRFP = imfinfo(RFPfilename);
for ii = 2:length(tiff_infoRFP)
    temp_tiffRFP = imread(RFPfilename, ii);
    tiff_stackRFP = cat(3 , tiff_stackRFP, temp_tiffRFP); % RFP inensity 
end


for znum = 1:size(tiff_stack,3)
    gfpPixelListeach = bfPixelList(bfnum,1);
    nucintensity = coor2inten(gfpPixelListeach,tiff_stack(:,:,znum));
    nucintensityGFP(znum,1:3) = nucintensity; % get the gfp BF seg intensity 

    nucintensityRFP = coor2inten(gfpPixelListeach,tiff_stackRFP(:,:,znum));
    nucintensity2RFP(znum,1:3) = nucintensityRFP; % get the RFP BF seg intensity 
end
%%
maxintenGFPdou = cell2mat(nucintensityGFP(:,3));
maxintenRFPdou =cell2mat(nucintensity2RFP(:,3));
%%
norintenGFP = (ra.* maxintenGFPdou - maxintenRFPdou)./(ra - rg);
norintenGFPall = mean(norintenGFP).*15;
%%
a = ischange(maxintenGFPdou,"linear","MaxNumChanges",3);
b = find(a == 1);
maxintenGFPdou2 = maxintenGFPdou(min(b):max(b),:);
maxintenRFPdou2 = maxintenRFPdou(min(b):max(b),:);
%%
norintenGFP2 = (ra.* maxintenGFPdou2 - maxintenRFPdou2)./(ra - rg);
norintenGFPall2 = sum(norintenGFP2);
end
%%
function cellinten2 = coor2inten(nulPL,mipGFP)
cellinten = cell(1,1);
for cellnumber = 1:length(nulPL)
    celleachintensity = cell(1,1);
    celleachcoor = nulPL{cellnumber,1};
    for celleachsize = 1:length(celleachcoor)
        celleachintensity{celleachsize,1} = mipGFP(celleachcoor(celleachsize,2),celleachcoor(celleachsize,1));
    end
    cellinten{cellnumber,1} = cell2mat(celleachintensity);
    cellsuminten{cellnumber,1} = sum(cell2mat(celleachintensity));
end
cellinten2 = horzcat(nulPL,cellinten,cellsuminten);
end


% GFP channel nuclear mask 
function [GFP,nulPL,tiff_stack] = mask_GFP_copynum2(GFPfilename,GFP_seg,savename)
tiff_stack = imread(GFPfilename, 1) ; % read in first image
%concatenate each successive tiff to tiff_stack
tiff_info = imfinfo(GFPfilename);

for ii = 2:length(tiff_info)
    temp_tiff = imread(GFPfilename, ii);
    tiff_stack = cat(3 , tiff_stack, temp_tiff);
end
GFP = GFP_seg;
stats = regionprops(GFP,'PixelList');
nulPL = {stats.PixelList}.';
end

function volume = volume_estimate(bfmask)
mask = bfmask;
labels = max(unique(bfmask));
if labels == 0
    volume = []; 
else
for l = 1:labels %for each cell
    curr_mask = mask;
    
    % __ volume measurements below __
    curr_mask(curr_mask ~= l) = 0; %leave only mask for current cell
    curr_mask = curr_mask ./ l; % convert mask value to 1
    curr_mask_RP=regionprops(curr_mask,'boundingbox','Orientation'); %get bounding box and major and minor axes
    
    xc=round(curr_mask_RP.BoundingBox(1)):-1+round((curr_mask_RP.BoundingBox(1)+curr_mask_RP.BoundingBox(3))); %crop mask to the bounding box
    yc=round(curr_mask_RP.BoundingBox(2)):-1+round((curr_mask_RP.BoundingBox(2)+curr_mask_RP.BoundingBox(4)));
    curr_mask1=curr_mask(yc,xc);
    
    rotated_mask = imrotate(curr_mask1,-curr_mask_RP.Orientation); %rotate cell mask so that its major axis is horizontal
    
    rotated_mask_RP=regionprops(rotated_mask,'boundingbox'); %crop mask to bounding box again
    xc2=round(rotated_mask_RP.BoundingBox(1)):-1+round((rotated_mask_RP.BoundingBox(1)+rotated_mask_RP.BoundingBox(3)));
    yc2=round(rotated_mask_RP.BoundingBox(2)):-1+round((rotated_mask_RP.BoundingBox(2)+rotated_mask_RP.BoundingBox(4)));           
    rotated_mask=rotated_mask(yc2,xc2);
    
    [xs,ys]=size(rotated_mask);
    vol_ht = NaN(1,ys); 
    for i=1:ys; %for each pixel (aka unit distance) in major axis
        r_h=(sum(rotated_mask(:,i))./2); %find the mean distance of both orthogonals above and below the pixel to the edge of the mask
        vol_ht(i)=pi.*(r_h).^2; %use that value as the radius to find the area of a circle
    end
    volume(l) = sum(vol_ht); % sum all the areas together to get an estimate of the volume
end
end
end

function [bf_label, bfPixelList] = mask_manual_BF_auto_seg(BF,savename)
 BF = ~imbinarize(BF);
BF_file4 = imfill(BF,'holes');
bf_label = bwlabel(BF_file4,8);
stats = regionprops(bf_label,'PixelList');
bfPixelList = {stats.PixelList}.';

save_cell = strrep(savename,'BF','bfPixelList');
save(save_cell, 'bfPixelList');
save_image = strrep(savename,'BF','manul_mask');
save(save_image, 'bf_label');
end
