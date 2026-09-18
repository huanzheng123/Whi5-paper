% calculate the diffusion coefficient from track data 
% get the trackinfo first, combine all the useful tracks in a same file,
% then use specific MSD method calculate the D
function [Dcell,trackcell,D_alpha_sig] = cal_D_1stMSD(trackdata,DhistMinSteps,savename,D1,loci,locl,locu)
% trackfile2 = uigetfile('*tracks_record*.mat', 'Multiselect', 'on','Pick your filtered track data');
% trackdata = trackfile2;
%%
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
% Apply MSD to the track file 
datasize = size(trackinfo_sum_fil);
D_alpha_sig = cell(datasize);
D_alpha_sig2 = cell(datasize);
for i = 1:datasize(1)
        for j = 1:datasize(2)
            celltracks = trackinfo_sum_fil{i,j};
            D_alpha_sig_cell = cell(size(celltracks));
           
            for k = 1:length(celltracks)
                tracks = celltracks{k,1};
                [para] = MSD_cal(tracks,2);
                D_alpha_sig_cell{k} = para;
                
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
trackcell = trackinfo_sum_fil;

