dir_name = "AWA015_PTA_1_Rec_Trans/downsampled/muscle_segmentation/";
extension = "png";
regions = ["left", "right"];
nb_used_slices = 5;
base_dir = join([getenv("HOME"), "Documents/phd/microCT/data", dir_name], '/');

toml_map = toml.read(base_dir + "analysis.toml");
params = toml.map_to_struct(toml_map);

mask_paths = getImagePaths(base_dir, extension);
mask_stack = loadImageStack(mask_paths);

for k = 1:length(regions)
    region = regions(k);

    if strcmp(region, "left")
        start_nb = params.left.start_nb;
        end_nb = params.left.end_nb;
    elseif strcmp(region, "right")
        start_nb = params.right.start_nb;
        end_nb = params.right.end_nb;
    elseif strcmp(region, "body")
        start_nb = params.body.start_nb;
        end_nb = params.body.end_nb;        
    else
        error("Error: invalid horn selection.");
    end

    disp("Rotating region: " + region);
    rotated_stack = rotateImageStack( ...
        mask_stack(:, :, start_nb:end_nb), region, nb_used_slices); 

    disp("Saving region: " + region);
    saveImageStack(rotated_stack, base_dir + region, ...
        params.prefix, 0, extension);

    clear rotated_stack % Save memory
end
