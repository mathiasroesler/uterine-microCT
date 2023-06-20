function uCTCentreline(dir_name, regions, extension)
%UCTCENTRELINE Computes the centreline of the images found in the provided
%directory based on the selected regions.
%
%   The base directory is $HOME/Documents/phd
%
%   Input:
%    - dir_name, path to the dataset from the base directory, should end 
%    with a /.
%    - regions, name of the regions that are going to be processed, either
%    left, right of [left, right], default value is [left, right].
%    - extension, extension of the image, default value is png
%
%   Return:
if nargin < 3
    extension = "png";
end

if nargin < 2
    regions = ['left', 'right'];
end

base_dir = join([getenv("HOME"), "Documents/phd/", dir_name], '/');

for k = 1:length(regions)
    region = regions(k);
    disp("Processing region: " + region)
    mask_paths = getImagePaths(base_dir + region, extension);
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

    save(base_dir + region + "/centreline.mat", "centreline");
    clear centreline
end
end