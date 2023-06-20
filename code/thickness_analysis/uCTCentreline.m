function uCTCentreline(dir_path, base_name, regions, downsampled, extension)
%UCTCENTRELINE Computes the centreline for the dataset provided by
%base_name given the selected regions. The funciton assumes that the masks
%that are used to compute the centreline are located in a
%muscle_segmentation folder.
%   
%   base_dir is $HOME/Documents/phd/
%
%   Input:
%    - dir_path, path to the directory containing the dataset from base_dir
%    - base_name, name of the dataset.
%    - region, either left, right or both, used to sort the centrepoints. 
%    - downsampled, true if the dataset has been downsampled, default value
%    is true.
%    - extension, extension of the images to load, default value is png.
%
%   Return: 
if nargin < 5
    extension = "png";
end
if nargin < 4
    downsampled = true;
end

base_dir = join([getenv("HOME"), "Documents/phd/", dir_path, base_name], '/');

if downsampled
    base_dir = join([base_dir, "downsampled"], '/');
end

base_dir = join([base_dir, "muscle_segmentation"], '/');

for k = 1:length(regions)
    region = regions(k);
    disp("Processing region: " + region)

    if strcmp(region, "both")
        mask_paths = getImagePaths(base_dir, extension);
    else
        mask_paths = getImagePaths(base_dir + region, extension);
    end

    mask_stack = loadImageStack(mask_paths);

    nb_slices = size(mask_stack, 3);
    centreline = zeros(6, nb_slices); % Placeholder for 3 centre points

    for m = 1:nb_slices
        centre_points = findCentrepoints(mask_stack(:, :, m), region);

        if size(centre_points, 1) == 3
            centreline(:, m) = reshape(centre_points', [6, 1]);
        else
            if matches(region, "left")
                centreline(1:2, m) = centre_points;
            elseif matches(region, "right")
                centreline(1:2, m) = centre_points;
            end
        end
    end

    if strcmp(region, "both")
        save(base_dir + "/centreline.mat", "centreline");
    else
        save(base_dir + region + "/centreline.mat", "centreline");
    end
    clear centreline
end
