function removeOvaries(dir_path, base_name, extension, start_nb)
%REMOVEOVARIES Removes the ovaries from the segmentation masks of a uCT
%dataset.
%
%   base_dir is $HOME/Documents/phd/
%
%   Input:
%    - dir_path, path to the directory containing the dataset from base_dir
%    - base_name, name of the dataset.
%    - extension, extension of the segmentation masks, default png.
%    - start_nb, number at which to start saving images, default 1.
%   Return:
if nargin < 5
    start_nb = 1;
end

if nargin < 4
    extension = "png";
end


load_directory = join([getenv("HOME"), "Documents/phd", dir_path, base_name], '/');
img_paths = getImagePaths(load_directory, extension);

% Load parameters
toml_map = toml.read(join([load_directory, base_name + ".toml"], '/'));
params = toml.map_to_struct(toml_map);

disp('Loading image stack')
img_stack = loadImageStack(img_paths);

middle_pixel = floor(max(size(img_stack)) / 2);
[~, largest_dim] = max(size(img_stack));
nb_img = size(img_stack, 3);

disp('Removing left ovary')
for k = params.downsampled.left_ovary:nb_img
    % Remove ovary on the left side of the image
    switch largest_dim
        case 1
            img_stack(1:middle_pixel, :, k) = 0;

        case 2
            img_stack(:, 1:middle_pixel, k) = 0;
    end
end

disp('Removing right ovary')
for k = params.downsampled.right_ovary:nb_img
    % Remove ovary on the right side of the image
    switch largest_dim
        case 1
            img_stack(middle_pixel:end, :, k) = 0;

        case 2
            img_stack(:, middle_pixel:end, k) = 0;
    end
end

disp('Saving image stack')
img_prefix = loadParams(load_directory + '/' + params_file);
saveImageStack(img_stack, load_directory, img_prefix, start_nb, extension);

end