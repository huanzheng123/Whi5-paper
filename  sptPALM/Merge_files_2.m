% merge BF and GFP file 
clear;close all;
BFfiles = "dwell_45_20220223_MS162_200ms_BF_0%d.tif";
GFPfiles = "dwell_45_20220223_MS162_200ms_GFP_0%d.tif";
fieldview = [1:1:7];

% read files 
BFlist = readfiles(fieldview,BFfiles);
GFPlist = readfiles(fieldview, GFPfiles);

% merge files 
% GFPmerge(BFlist,GFPlist)



%function GFPmerge(BFlist,GFPlist)
for j = 1:length(BFlist)
    BFfile = imread(BFlist{j});
    GFPimage = imread(GFPlist{j},1);
    tiff_stack = [];
    tiff_info = imfinfo(GFPlist{j}); % return tiff structure, one element per image
for ii = 2 : size(tiff_info, 1)
        temp_tiff = imread(GFPlist{j}, ii);
    tiff_stack = cat(3 , tiff_stack, temp_tiff);
end
GFPmip = max(tiff_stack, [], 3);
empty = zeros(size(GFPmip,1),size(GFPmip,2));
% 
% BFimage2 = (BFfile-min(BFfile))./(max(BFfile) - min(BFfile));
% GFP_MAX2 = (GFPmip-min(GFPmip))./(max(GFPmip) - min(GFPmip));

GFP_MAX2 = mean(mean(BFfile))./mean(mean(GFPmip)).*GFPmip;
GFP_MAX22 = mat2gray(GFPmip);
BFfile2 = mat2gray(BFfile);

newmerge = cat(3,GFP_MAX22,BFfile2,BFfile2);
    figure
    image(newmerge);
   savename = strrep(GFPlist{j},'GFP','merge_GFP');
   savename2 = strrep(savename,'.tif','.png');
   imwrite(newmerge,savename2);
    
end

% end
%%
function filelist = readfiles(fieldview, filebasename)
filelist = cell(1,1);
for i = 1:length(fieldview)
filename = sprintf(filebasename,fieldview(i));
filelist{i} = filename;
end
end
