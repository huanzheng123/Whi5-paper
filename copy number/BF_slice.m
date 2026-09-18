% scipt the substack the in focus channel 
% slice the BF z stack 
fileall = [2 3];
slicenum = [12];

for fileiter = 1:length(fileall)
if fileall(fileiter)<10
    BF_filename = sprintf('97_20240116_YHZ60_0%d_w1BF.tif',fileall(fileiter));
else
    BF_filename = sprintf('20230923_FISH_ZEY281_01_%d_BF.tif',fileall(fileiter));
    
end
BF_file= imread(BF_filename,slicenum);
savename = strrep(BF_filename,'BF','_sliced_BF');
imwrite(BF_file,savename);

end

