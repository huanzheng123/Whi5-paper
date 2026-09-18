% Specify the folder containing the binary images
clear;close all; 
iter = [11];
for iii = 1:length(iter)
imageFolder = '/Users/zhenghuan/Downloads/files for Masoumeh/YHZ75/66';
filename = sprintf('dwell_66_20230915_YHZ75_%d_w2T2_merge_manul_mask_gfp-*',iter(iii));
imageFiles = dir(filename); % Assuming .mat files

savename = strrep(imageFiles(1).name,'merge_manul_mask_gfp','');
savename2 = strrep(savename,'mat','png');
% Initialize an empty array for the combined image
combinedImage = zeros(512, 512); % Assuming all images are 512x512

% Loop through each file, load it, and combine
for i = 1:length(imageFiles)
    % Load the current image
    filePath = fullfile(imageFiles(i).name);
    data = load(filePath);
    
    % Assuming the binary image variable is named 'binaryImage' in the .mat file
    fieldNames = fieldnames(data);
    binaryImage = data.(fieldNames{1}); % Retrieve the image
    
    % Combine: You can use + for sum or | for logical OR
    combinedImage = combinedImage + binaryImage; % For sum
    % combinedImage = combinedImage | binaryImage; % For logical OR
end

% Normalize if summing results in values >1 (optional)
% combinedImage = combinedImage / max(combinedImage(:));

% Display the result
outputFileName = fullfile(imageFolder, savename2);
imwrite(combinedImage, outputFileName); % Save as PNG

end
