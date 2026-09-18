% GFP intensity extraction 
clear, close all;
ra = 0.8798; % use this ra for now, need to estimate ra from MS10 later; 
rg = 0.1065; %from pure mNeonGreen data 
ratiov = 0.0995; % ratio of cell volume and nuclear volume in BY4741; 
pixelsize = 0.0725; % um per pixel; 

% input data: BF mask and GFP z 
% load the different BF mask layer
fileall = [0:1:14];
cellpro = [3 4 5 6 7]; 

for fileiter = 1:length(fileall)
    for cellproiter = 1:length(cellpro)
        close all;
        fileidx = fileall(fileiter);
        celltypeiter = cellpro(cellproiter);

if fileidx < 100
    BF_filename = sprintf('20231218_PAY09_%d_w2merge_GFP-%d.mat',fileidx,celltypeiter);
    BF_file= imread(BF_filename);
    BF_file= BF_file(:,:,1);
    BF_oriname = sprintf('20231218_PAY09_%d_w2_BF_Cropped.TIF',fileidx);
    BF_ori = imread(BF_oriname);
    savename = strrep(BF_filename,'GFP','BF');
    GFPfilename = sprintf('20231218_PAY09_%d_w1_Cropped_GFP-488-561-610-75new.TIF',fileidx);
    RFPfilename = strrep(GFPfilename,'GFP','RFP');
else
    BF_filename = sprintf('20220717_MS564_%d_w1merge_GFP-%d.mat',fileidx,celltypeiter);
    BF_file= imread(BF_filename);
    BF_file= BF_file(:,:,1);
    BF_oriname = sprintf('20220717_MS564_%d_w1_BF_Cropped.TIF',fileidx);
    BF_ori = imread(BF_oriname);
    savename = strrep(BF_filename,'GFP','BF');
    GFPfilename = sprintf('20220717_MS564_%d_w2_Cropped_GFP-488-561-610-75.TIF',fileidx);
    RFPfilename = strrep(GFPfilename,'GFP','RFP');
end

    


  [BFManuMask,bfPixelList] = mask_manual_BF_newSeg(BF_file,savename);
  [GFPManuMask,gfpPixelList,tiff_stack] = mask_GFP_copynum(GFPfilename,savename);

  volume = volume_estimate(BFManuMask);
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
areaind = find(Darea2(:,3)==0);
Darea2(areaind,:)=[];

% get the slice number and intensity for inquired GFP cell
for GFPiter = 1:size(Darea2,1)
    GFPnum = Darea2(GFPiter,2);
    [slicenum, maxintenGFP,maxintenRFP] = zlayerinten(tiff_stack, RFPfilename, gfpPixelList,GFPnum); 
    Darea2(GFPiter,6) = slicenum;
    Darea2(GFPiter,7) = maxintenGFP;
    norintenGFP = (ra.* maxintenGFP - maxintenRFP)./(ra - rg);
    Darea2(GFPiter,8) = norintenGFP; % integreted intensity of the slice of nuclear 
    wholenucinten = norintenGFP./Darea2(GFPiter,5).*0.0995.* Darea2(GFPiter,4); 
    Darea2(GFPiter,9) = wholenucinten;
    Darea2(GFPiter,10) = (wholenucinten+201.96)./1742.1; % y = 1742.1X-201.96 mNG copy num relation with intensity;
    Darea2(GFPiter,11) = Darea2(GFPiter,10)./Darea2(GFPiter,4);
end

% save cell area: BF labeled mask and GFP labeled mask; 1.BF idx 2.nuclear idx 3. BF
% cell area 4. BF cell v 5. nuclear area 6.max intensity slice num 7. intensity of
% that slice 8.normailized intensity 9.estimated whole v intensity; 10.
% normalized copy number 11. concentration 

savename2 = strrep(savename,'BF','intensity_GFP');
savename3 = savename2()
save(savename2,'Darea2');


    end
end

% structure of the data, BF idx, BF mask, BF cell size and volume, GFP cell
% idx, nuclear size and intensity, normalized intensity 

%%
% input coordinates and image, output intensity 
function [slicenum, maxintenGFP,RFPinten] = zlayerinten(tiff_stack, RFPfilename, gfpPixelList,GFPnum)
% read RFP channel stack 
tiff_stackRFP = imread(RFPfilename, 1) ; % read in first image
%concatenate each successive tiff to tiff_stack
tiff_infoRFP = imfinfo(RFPfilename);
for ii = 2:length(tiff_infoRFP)
    temp_tiffRFP = imread(RFPfilename, ii);
    tiff_stackRFP = cat(3 , tiff_stackRFP, temp_tiffRFP);
end



for znum = 1:size(tiff_stack,3)
    gfpPixelListeach = gfpPixelList(GFPnum,1);
    nucintensity = coor2inten(gfpPixelListeach,tiff_stack(:,:,znum));
    nucintensity2(znum,1:3) = nucintensity;

    nucintensityRFP = coor2inten(gfpPixelListeach,tiff_stackRFP(:,:,znum));
    nucintensity2RFP(znum,1:3) = nucintensityRFP;
end
% find the slice that have the max intensity for each nuclear
[maxintenGFP,slicenum] = max(cell2mat(nucintensity2(:,3)),[],'all');
% find the intensity of slicenum 
RFPinten = nucintensity2RFP{slicenum,3};
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
function [GFP,nulPL,tiff_stack] = mask_GFP_copynum(GFPfilename,savename)
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
image_bi = imbinarize(bw,'adaptive',Sensitivity=0.2);
imagecell = bwareaopen(image_bi,200);
image_label = bwlabel(imagecell,8);
imshowpair(mipGFP,imagecell);

stats = regionprops(image_label,mipGFP,'Image','PixelList','PixelValues');
PixelValues = {stats.PixelValues}.';
image_binary = cell(length(PixelValues),1);
PixelValues = {stats.PixelValues}.';
PixelList = {stats.PixelList}.';
Pixel_binary = cell(length(PixelList),1);
% save GFP_mask_01_filtered.mat image_label;
GFP = image_label;
nulPL = PixelList;
% save(savename,"GFP");
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
