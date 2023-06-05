dir_name = "AWA015_PTA_1_Rec_Trans/downsampled/muscle_segmentation/";
extension = "png";
horns = ["left", "right"];
base_dir = join([getenv("HOME"), "Documents/phd/microCT/data", dir_name], '/');

for k = 1:length(horns)
    horn = horns(k);
    mask_paths = getImagePaths(base_dir + horn + "_horn", extension);
    mask_stack = loadImageStack(mask_paths);

    nb_slices = size(mask_stack, 3);
    centreline = zeros(2, nb_slices);

    for m = 1:nb_slices
        centre_points = findCentrepoints(mask_stack(:, :, m));
        centreline(:, m) = centre_points(k, :);
    end

    save(base_dir + horn + "_horn/centreline.mat", "centreline");
    clear centreline
end
