files = dir('*.png');
for ii = 1:length(files)
    % Get the file name (minus the extension)
    [~, fname] = fileparts(files(ii).name);
    movefile(files(ii).name, sprintf('%s.mat',fname));
end