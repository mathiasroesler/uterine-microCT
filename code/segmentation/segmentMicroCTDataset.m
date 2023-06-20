function segmentMicroCTDataset(dir_path, base_name, segmentation_type, ...
    split_nb, downsampled, extension, start_nb)
%SEGMENTMICROCTDATASET Segments a uCT dataset based on the segmentation 
%type. 
%
%   base_dir is $HOME/Documents/phd/
%
%   Input:
%    - dir_path, path to the directory containing the dataset from base_dir
%    - base_name, name of the dataset.
%    - segmentation_type, type of segmentation to use:
%      1 is tissue segmentation
%      2 is muscle layer segmentation
%      3 is fat segmentation
%      4 is shape segmentation
%    - split_nb, slice at which to split into left and right horn.
%    - downsampled, true if the dataset has been downsampled, default value
%    is true.
%    - extension, extension of the segmentation masks, default png.
%    - start_nb, number at which to start saving images, default 1.
%   Return:
if nargin < 7
    start_nb = 1;
end
if nargin < 6
    extension = "png";
end
if nargin < 5
    downsampled = true;
end

load_directory = join([getenv("HOME"), "Documents/phd", dir_path, base_name], ...
    '/'); % Directory where images are located

% Load parameters
toml_map = toml.read(join([load_directory, base_name + ".toml"], '/'));
params = toml.map_to_struct(toml_map);
preprocess = params.preprocess;
morph_size = params.morph_size;
sigma = params.sigma;

if downsampled
    % If using the downsampled dataset
    load_directory = join([load_directory, "downsampled"], '/');
    preprocess = params.downsampled.preprocess;
    morph_size = params.downsampled.morph_size;
    sigma = params.downsampled.sigma;
end

img_paths = getImagePaths(load_directory, extension);

disp('Loading image stack')
img_stack = loadImageStack(img_paths);

disp('Spliting image stack')
stack_cell = splitImageStack(img_stack, split_nb);
mask_cell = cell(size(stack_cell)); % Empty cell array for the masks

for k = 1:length(stack_cell)
    img_stack = stack_cell{k};
    disp('Processing stack ' + string(k))
    if preprocess ~= 0
        disp('    Preprocessing image stack')
        if preprocess == 4
            mask_paths = getImagePaths(load_directory + "/tissue_segmentation", ...
                extenstion);
            mask_stack = loadImageStack(mask_paths);
            img_stack = preprocessImageStack(img_stack, preprocess, sigma, ...
                mask_stack);
            clear mask_paths mask_stack;
        else
            img_stack = preprocessImageStack(img_stack, preprocess, ...
                morph_size);
        end
    end

    disp('    Segmenting image stack')
    switch segmentation_type
        case 1
            % Tissue segmentation
            save_dir= join([load_directory, "tissue_segmentation"], '/');
            mask_stack = tissueSegmentation(img_stack, ...
                params.segmentation.tissue.threshold, ...
                params.segmentation.tissue.morph_size, ...
                params.segmentation.tissue.neighborhood_size);

        case 2
            % Muscle segmentation
            save_dir = join([load_directory, "muscle_segmentation"], '/');
            mask_stack = muscleSegmentation(img_stack, ...
                params.segmentation.muscle.morph_size, ...
                params.segmentation.muscle.neighborhood_size);

        case 3
            % Fat segmentation
            save_dir= join([load_directory, "fat_segmentation"], '/');
            mask_stack = fatSegmentation(img_stack, ...
                params.segmentation.fat.morph_size);

        case 4
            % Shape segmentation
            save_dir= join([load_directory, "shape_segmentation"], '/');
            mask_stack = shapeSegmentation(img_stack, ...
                params.segmentation.shape.threshold, ...
                params.segmentation.shape.morph_size, ...
                params.segmentation.shape.neighborhood_size);
    end

    mask_cell{k} = mask_stack;
end

disp('Fusing mask stack')
mask_stack = fuseImageStacks(mask_cell);

disp('Saving mask stack')
saveImageStack(mask_stack, save_dir, params.img_prefix, start_nb, extension);


end