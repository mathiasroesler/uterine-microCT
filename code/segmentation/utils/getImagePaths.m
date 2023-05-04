function img_paths = getImagePaths(dir_path, extension)
%GETIMAGESPATHS Retrieves the paths to the images in the provided folder
%with the given extension.
%   
%   Input:
%    - dir_path, path to folder containing the images.
%    - extension, file extension to search for in folder, default .bmp.
%   
%   Return:
%    - img_paths, cell of paths to the images.
if nargin < 2
    extension = ".bmp";
end

char_extension = char(extension);

if char_extension(1) ~= '.'
    extension = "." + extension;
end

listings = dir(dir_path);
len = length(listings);

img_paths = cell(len, 1);
home_dir = convertCharsToStrings(listings(1).folder);

for k = 1:len % Retrieve only file with the correct extension
    [~, ~, file_extension] = fileparts(listings(k).name);

    if ~listings(k).isdir && file_extension == extension
        name = convertCharsToStrings(listings(k).name);
        img_paths{k} = join([home_dir, name], '/');
    end
end

img_paths = img_paths(~cellfun('isempty', img_paths)); % Remove empty cells

end