% get the dwell time estimate in single cells with Whi5 as nuclear marker 
% run this after the spt_Whi5_nuclear script 
clear; close all;
% new version, add the cell volume to the end of the Dareacell column 
DhistMinSteps = 4;
Rthre = 2;
Nthre = 4;
jumpgap = 1; % 1 step

fileall = [1 2 3 5]; % input the index of field of view 
celltypeall = [3:1:7];
for celltypeaiter = 1:length(celltypeall)
    celltypeiter2 = celltypeall(celltypeaiter);
for fileiter = 1:length(fileall)
fileidx = fileall(fileiter);
BF_filename = sprintf('dwell_68_20230918_YHZ68_0%d_w2T2_merge_GFP-%d.mat',fileidx,celltypeiter2);
BF_filename = strrep(BF_filename,'.tif','.mat');
savename = BF_filename;  
save_cellarea = strrep(BF_filename,'GFP','cell_area');
cellarea = load(save_cellarea);
cellarea = cellarea.cellarea;
bffile = strrep(BF_filename,'GFP','manul_mask');
gfpfile = strrep(BF_filename,'GFP','gfpPixelList');
gfpPixelList = load(gfpfile);
gfpPixelList = gfpPixelList.gfpPixelList;
% read track file 
trackfilename = sprintf('Stream%d_dwell_68_20230918_YHZ68_0%d_w3T3 RFP CUST.tif_spots.csv',fileidx,fileidx);
trackfiledir = dir(trackfilename);
trackfile = cell(1,1);
for trackfileiter = 1:length(trackfiledir)
temptrack = trackfiledir(trackfileiter);
trackfile{trackfileiter} = temptrack.name;
end

bfmask = load(bffile);
bfmask = bfmask.bf_label;
% get the cell volumn info 
if max(max(bfmask)) == 0
    cellvolume = [];
else
cellvolume = volume_estimate(bfmask);
end
%%
% link to the old script 
trackfile2 = track_seg_function(trackfile,gfpPixelList); % filter the tracks outside the nucleus area
[Dcell,trackcell,cellnum] = cal_D_1stMSD(trackfile2,DhistMinSteps,savename); % linked to the old file without chagning the data format; 
%returin the dwell time results, in a strucutre format; 
Dwellres = cal_dwell(trackfile2,DhistMinSteps,Nthre,Rthre,jumpgap,savename);

% add the index of trackcell into Dcell matrix 
DcellInd = Dcell;
DcellSize = size(Dcell);
for i = 1:DcellSize(1)
    for j = 1:DcellSize(2)
        DcellDou = Dcell{i,j};
        Ind1 = repmat(i,length(DcellDou),1);
        Ind2 = repmat(j,length(DcellDou),1);
        Ind3 = [1:1:length(DcellDou)]';
        DcellDou2 = [DcellDou Ind1 Ind2 Ind3];
        DcellInd{i,j} = DcellDou2;
    end
end

cellnum2 = size(cellnum);
Dcell2 = {};
for i = 1:cellnum2(1)
    Dcell2{i} = vertcat(DcellInd{i,:});
end

Darea = [num2cell(cellarea(:,1)),num2cell(cellarea(:,2)),num2cell(cellarea(:,3)), Dcell2'];

% Delete the zero cell area
Darea2 = Darea;
Dareacol = cell2mat(Darea2(:,3));
areaind = find(Dareacol==0);
Darea2(areaind,:)=[];
Dwellres(areaind,:)=[];

% add the field index 
fieldinx = fileidx;
Darea2 = [[repmat({fieldinx},1,size(Darea2,1))' Darea2]];
for Darea2iter = 1:size(Darea2)
    Darea2{Darea2iter,4} = [Darea2{Darea2iter,4},cellvolume(1,cell2mat(Darea2(Darea2iter,2)))];
end

Darea3 = horzcat(Darea2,Dwellres);
% add the field index to the track info file 
trackcell2 = [trackcell repmat({fieldinx},size(trackcell,1),1)];
% add the cell index into the track info file 

trackcellinx = [1:1:size((trackcell2),1)]';
trackcell3 = [trackcell2 num2cell(trackcellinx)];

% save the data 
RNAP2Result = struct();
settings = struct();
% settings.pixel = pixel; 
% settings.dT = dT;
% settings.sigmaNoise = sigmaNoise;
settings.DhistMinSteps = DhistMinSteps;
Result.settings = settings;
Result.Darea = Darea2;
save_Darea2 = strrep(savename,'GFP','Darea_Dwellres');
save(save_Darea2,'Darea3');
save_WS = strrep(savename,'GFP','WS_RG');
save(save_WS);
savetrack = strrep(savename,'GFP','trackinfo_RG');
save(savetrack,'trackcell3');
end
end
% function to convert BF mask to volumn 
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

function [Dwellres] = cal_dwell(trackdata,DhistMinSteps,Nthre,Rthre,jumpgap,savename)
num_files = length(trackdata);
Dcell = {};
trackinfo_1 = importdata(trackdata{1});
trackinfo_sum = cell(length(trackinfo_1),num_files);
%%
for z = 1:num_files %iterate each track mat file
    a = importdata(trackdata{z});
    for j = 1:length(a)
        b = a{j};
        trackinfo = b(~all(cellfun('isempty', b(:, 1)), 2), :);
        trackinfo_sum(j,z) = {trackinfo};
    end
end
%%
% filter out the short tracks in trackinfo file
trackinfo_sum_fil = trackinfo_sum;
trackinfosize = size(trackinfo_sum);
for t = 1:trackinfosize(1)
    for r = 1:trackinfosize(2)
        W = trackinfo_sum_fil{t,r};
        W = W(cellfun(@(x) length(x) >= DhistMinSteps, W));
        trackinfo_sum_fil{t,r} = W;
    end
end
%%
% Apply dwell time analysis to the track file 
datasize = size(trackinfo_sum_fil);
D_alpha_sig = cell(datasize);
D_alpha_sig2 = cell(datasize);
for i = 1:datasize(1)
        for j = 1:datasize(2)
            celltracks = trackinfo_sum_fil{i,j};
            D_alpha_sig_cell = cell(size(celltracks));
           
            for k = 1:length(celltracks)
                tracks = celltracks{k,1};
                trackresults = classify_tracks(tracks,Nthre,Rthre,jumpgap);
                D_alpha_sig_cell{k} = trackresults;
                
            end  
            D_alpha_sig{i,j} = D_alpha_sig_cell;
           
         end
end
%%
Dcell = cell(datasize);
for i = 1:datasize(1)
        for j = 1:datasize(2)
            c = D_alpha_sig{i,j};
            Dcell{i,j} = cellfun(@(v)v(1),c);
        end
end

Dwellres = Dcell;

end