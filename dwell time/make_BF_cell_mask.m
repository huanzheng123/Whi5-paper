% make the GFP mask file based on cell size 
% from the BF mask get the volume info -- from the volume set the
% threshold -- filter out the coupled GFP info; 

filename = '*merge_manul_mask*%d.mat';
cellsizethre = [20,50]; 

% input the GFP and BF mask files 
Dareafilename = sprintf('*merge_Darea_Dwellres*%d.mat',celltype);



% 
% bffile = strrep(BF_filename,'GFP','manul_mask');
% bfmask = load(bffile);
% bfmask = bfmask.bf_label;
bfmask = bf_label;
% get the cell volumn info 
if max(max(bfmask)) == 0
    cellvolume = [];
else
cellvolume = volume_estimate(bfmask);
end

cellvol = cellvolume.*0.10^3;
% apply the cellsize threshold 
a = find(cellvol < cellsizethre(1) & cellvol > cellsizethre(2));

if isempty(a) == 1
    aaa = 1;
else
    for i = 1:length(a)
        b = find(bfmask == a(i));
        bfmask(b) = 0;
    end
end
imshow(bfmask)













function volume = volume_estimate(bfmask)
mask = bfmask;
labels = max(unique(bfmask));
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