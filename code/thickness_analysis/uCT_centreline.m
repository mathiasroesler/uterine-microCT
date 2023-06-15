% dir_name = "histology/muscle_segmentation/";
dir_name = "microCT/data/AWA015_PTA_1_Rec_Trans/downsampled/muscle_segmentation/";
extension = "png";
regions = ["left", "right"];
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
