function uCTRotation(dir_path, base_name, regions, nb_used_slices, ...
    downsampled, extension)
%UCTCENTRELINE Computes the centreline for the dataset provided by
%base_name given the selected regions.
%   
%   base_dir is $HOME/Documents/phd/
%
%   Input:
%    - dir_path, path to the directory containing the dataset from base_dir
%    - base_name, name of the dataset.
%    - region, either left, right or both, used to sort the centrepoints. 
%    - nb_used_slices, number of slices to use to determine centre vector,
%    default value 5.
%    - downsampled, true if the dataset has been downsampled, default value
%    is true.
%    - extension, extension of the images to load, default value is png.
%
%   Return: 
if nargin < 6
    extension = "png";
end
if nargin < 5
    downsampled = true;
end
if nargin < 4
    nb_used_slices = 5;
end

load_directory = join([getenv("HOME"), "Documents/phd/", dir_path, base_name], '/');

% Load parameters
toml_map = toml.read(join([load_directory, base_name + ".toml"], '/'));
params = toml.map_to_struct(toml_map);
start_nb = params.thickness.start_nb;

if downsampled
    % Deal with downsampled dataset
    load_directory = join([load_directory, "downsampled"], '/');
end

mask_paths = getImagePaths(load_directory, extension);
mask_stack = loadImageStack(mask_paths);

for k = 1:length(regions)
    region = regions(k);

    if strcmp(region, "left")
        end_nb = params.thickness.left.end_nb;
    elseif strcmp(region, "right")
        end_nb = params.thickness.right.end_nb;

        error("Error: invalid horn selection.");
    end

    disp("Rotating region: " + region);
    rotated_stack = rotateImageStack( ...
        mask_stack(:, :, start_nb:end_nb), region, nb_used_slices); 

    disp("Saving region: " + region);
    saveImageStack(rotated_stack, load_directory + region, ...
        params.prefix, start_nb, extension);

    clear rotated_stack % Save memory
end
