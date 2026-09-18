% crop image, output BF focuse slice, BF gfp composed image and splide gfp
% and rfp Z stack 
% before run this script, input the RFP and GFP images to fiji to substract the background 
% filename 
clear;close all;
%BFslicekeep = 10; 
%BFfiles = "20230630_copynum_whi5_mNG_0%d_w2BF.TIF";
GFPfiles = "20230630_copynum_swi6_mNG_%d_w1dual-488-561-610-75new.tif";
fieldview = [10:1:46];
%BFlist = readfiles(fieldview,BFfiles);
GFPlist = readfiles(fieldview, GFPfiles);

%%
% convert the BF z to BF slice, save the left one 
%BFimages = cropBFslice(BFlist,BFslicekeep);
GFPimages = cropGFPZ(GFPlist);

% for j = 1:length(BFlist)
%     %BFfile = BFimages{j,1};
%     GFPmip = GFPimages{j,1};
%     empty = zeros(size(GFPmip,1),size(GFPmip,2));
%     GFP_MAX2 = mean(mean(BFfile))./mean(mean(GFPmip)).*GFPmip;
%     GFP_MAX22 = mat2gray(GFPmip);
%     BFfile2 = mat2gray(BFfile);
%     newmerge = cat(3,GFP_MAX22,BFfile2,BFfile2);
%     figure
%     image(newmerge);
%    savename = strrep(BFlist{j},'BF','merge_GFP');
%    savename2 = strrep(savename,'.TIF','.png');
%    imwrite(newmerge,savename2);
% end

close all;
% separate the RFP and GFP channel, save both 
function GFPimages = cropGFPZ(GFPlist)
GFPimages = cell(1,1);
for i = 1:length(GFPlist)
    Z_file = GFPlist{i};
    Z_stack = length(imfinfo(Z_file));
    savename = strrep(Z_file,'dual','Cropped');
    savename_RFP = strrep(savename,'Cropped','_Cropped_RFP');
    savename_GFP = strrep(savename,'Cropped','_Cropped_GFP');
    tiff_stack = [];
    for j = 1:Z_stack
        each_plane = imread(Z_file,j);
        [RFP_p,GFP_p] = Splitting(each_plane);
        imwrite(RFP_p,savename_RFP,'TIFF','WriteMode','append')
        imwrite(GFP_p,savename_GFP,'TIFF','WriteMode','append')

        tiff_stack = cat(3 , tiff_stack, GFP_p);
    end
    GFPimages{i,1} = max(tiff_stack, [], 3);

end
end



function BFimages = cropBFslice(BFlist,BFslicekeep)
BFimages = cell(1,1);
for i = 1:length(BFlist)
    Z_file = BFlist{i};
    Z_stack = length(imfinfo(Z_file));
    savename = strrep(Z_file,'BF','_BF_Cropped');
    %savename_RFP = strrep(savename,'Cropped','Cropped_RFP');
    %savename_GFP = strrep(savename,'Cropped','Cropped_GFP');

        each_plane = imread(Z_file,BFslicekeep);
        [RFP_p,GFP_p] = Splitting(each_plane);
        imwrite(RFP_p,savename,'TIFF')
        BFimages{i,1} = RFP_p;
end
end


function filelist = readfiles(fieldview, filebasename)
filelist = cell(1,1);
for i = 1:length(fieldview)
filename = sprintf(filebasename,fieldview(i));
filelist{i} = filename;
end
end

function [RFPIma, GFPIma] = Splitting(orimage)
% ori_file = uigetfile('*Dual_Z*.tif', 'Multiselect', 'on');
% for i = 1:length(ori_Z)
% image = imread(orimage);
[~, columns, ~] = size(orimage);
midcol = ceil(columns/2);
GFPIma = orimage(:,midcol+1:end);
RFPIma = orimage(:,1:midcol);
end