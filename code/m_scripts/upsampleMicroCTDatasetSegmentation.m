function upsampleMicroCTDatasetSegmentation(dir_path, base_name, type, ...
    new_resolution, varargin)
%UPSAMPLEMICROCTDATASET Upsamples the segmentation of a downsampled uCT
%dataset. The new resolution is automatically calculated based on the
%config file of the non-downsampled dataset.
%
%   base_dir is $HOME/Documents/phd/ and set in utils/baseDir()
%
%   Input:
%    - dir_path, path to the directory containing the dataset from base_dir
%    - base_name, name of the dataset.
%    - type, segmentation type, {fat, tissue, shape, muscle}, default value
%    is muscle.
%    - new_resolution, either an array with the X, Y, and Z resolutions or
%    the factor by which to upsample by.
%    - batch_size, number of image to process in one batch, default 128.
%    - img_extension, extension of the images in the uCT dataset, default
%    png.
%   - save_extension, extension of the downsampled images, default png.
%
%   Return:
% Author Mathias Roesler July 2023
narginchk(4, 7);

if nargin < 7
    save_extension = "png";

else
    save_extension = varargin{3};
end

if nargin < 6
    img_extension = "png";

else
    img_extension = varargin{2};
end

if nargin < 5
    batch_size = 128;

else
    batch_size = varargin{1};
end

% Directory where images are located
main_directory = join([baseDir(), dir_path, base_name], '/');
img_load_directory = join([main_directory, "downsampled"], '/');
mask_load_directory = join([main_directory, "downsampled", ...
    type + "_segmentation"], '/');
img_save_directory = join([main_directory, "ST/masked"], '/');
mask_save_directory = join([main_directory, "ST/mask"], '/');

%% Load and set parameters
% Load properties of the original images that dont change
toml_map = toml.read(join([main_directory, base_name + ".toml"], '/'));
params = toml.map_to_struct(toml_map);

img_paths = getImagePaths(img_load_directory, img_extension);
mask_paths = getImagePaths(mask_load_directory, img_extension);
stack_size = length(img_paths);

if sum(size(new_resolution)) == 4
    % The input is the new resolution vector with 3 components
    new_nb_pixel_x = new_resolution(1);
    new_nb_pixel_y = new_resolution(2);
    new_stack_size = new_resolution(3);
    resize_factor = round(new_stack_size / stack_size);

elseif isscalar(new_resolution)
    % The input is the upsampling factory
    resize_factor = new_resolution(1);

else
    error("The size of the input resolution is wrong. Size should be 1 or 3.")
end

nb_runs = ceil(stack_size/batch_size); % Number of times to run loop
img_save_index = 1;

%% Main resizing loop
for run = 1:nb_runs
    disp("Running batch number " + num2str(run) + "/" + num2str(nb_runs));
    first_image_nb = (run-1) * batch_size + 1;
    last_image_nb = run * batch_size;

    if last_image_nb > stack_size
        last_image_nb = stack_size;
    end

    batch_stack_size = last_image_nb - first_image_nb + 1;
    new_batch_stack_size = round(batch_stack_size * resize_factor);

    % Load all images in current batch
    disp("Loading " + num2str(batch_stack_size) + " images in batch");
    img_stack = loadImageStack(img_paths(first_image_nb:last_image_nb));
    mask_stack = loadImageStack(mask_paths(first_image_nb:last_image_nb));

    % Resize current batch stack
    if sum(size(new_resolution)) == 4
        new_mask_stack = imresize3(uint8(mask_stack), [new_nb_pixel_x, ...
            new_nb_pixel_y, new_batch_stack_size]);
        new_img_stack = imresize3(uint8(img_stack), [new_nb_pixel_x, ...
            new_nb_pixel_y, new_batch_stack_size]);
    else
        new_mask_stack = imresize3(uint8(mask_stack), resize_factor);        
        new_img_stack = imresize3(uint8(img_stack), resize_factor);
    end
    
    % Apply masks
    new_img_stack = new_img_stack .* new_mask_stack;

    % Save current batch stack
    disp("Saving " + num2str(new_batch_stack_size) + " upsampled images");
    saveImageStack(new_img_stack, img_save_directory, params.prefix, ...
        img_save_index, save_extension);
    saveImageStack(imbinarize(new_mask_stack), mask_save_directory, ...
        params.prefix, img_save_index, save_extension);
    img_save_index = img_save_index + new_batch_stack_size;
end

log_file = join([img_load_directory, params.prefix + "_upsampled.log"], '/');
new_resolution = resize_factor*params.resolution; % um/pixel

% Log required information
file_ID = fopen(log_file, 'w');
fprintf(file_ID, "nb_pixel_x=%3d\n", size(new_img_stack, 1));
fprintf(file_ID, "nb_pixel_y=%3d\n", size(new_img_stack, 2));
fprintf(file_ID, "stack_size=%3d\n", stack_size*resize_factor);
fprintf(file_ID, "pixel_x_res=%03f\n", new_resolution);
fprintf(file_ID, "pixel_y_res=%03f\n", new_resolution);
fprintf(file_ID, "pixel_z_res=%03f\n", new_resolution);
% fclose(file_ID);