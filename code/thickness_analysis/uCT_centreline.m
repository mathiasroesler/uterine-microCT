dir_name = "AWA015_PTA_1_Rec_Trans/downsampled/muscle_segmentation/";
extension = "png";
regions = ["left", "right", "body"];
base_dir = join([getenv("HOME"), "Documents/phd/microCT/data", dir_name], '/');

for k = 1:length(regions)
    region = regions(k);
    mask_paths = getImagePaths(base_dir + region, extension);
    mask_stack = loadImageStack(mask_paths);

    nb_slices = size(mask_stack, 3);
    centreline = zeros(2, nb_slices);

    for m = 1:nb_slices
        centre_points = findCentrepoints(mask_stack(:, :, m), region);
        centreline(:, m) = centre_points(1, :);
    end

    save(base_dir + region + "/centreline.mat", "centreline");
    clear centreline
end
