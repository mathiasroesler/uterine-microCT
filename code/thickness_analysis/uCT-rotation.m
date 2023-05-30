dir_name = "AWA015_PTA_1_Rec_Trans/downsampled/muscle_segmentation/";
extension = "png";
horns = ["left", "right"];
nb_used_slices = 3;
base_dir = join([getenv("HOME"), "Documents/phd/microCT/data", dir_name], '/');

toml_map = toml.read(base_dir + "analysis.toml");
params = toml.map_to_struct(toml_map);

mask_paths = getImagePaths(base_dir, extension);
mask_stack = loadImageStack(mask_paths);

for k = 1:length(horns)
    horn = horns(k);

    if strcmp(horn, "left")
        start_nb = params.left.start_nb;
        end_nb = params.left.end_nb;
    elseif strcmp(horn, "right")
        start_nb = params.right.start_nb;
        end_nb = params.right.end_nb;
    else
        error("Error: invalid horn selection.");
    end

    rotated_stack = rotateImageStack( ...
        mask_stack(:, :, start_nb:end_nb), horn, nb_used_slices); 
    saveImageStack(rotated_stack, base_dir + horn + " _horn", ...
        params.prefix, 0, extension);

    clear rotated_stack % Save memory
end
