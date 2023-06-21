function downsampleMicroCTDataset(dir_path, base_name, varargin)
%DOWNSAMPLEMICROCTDATASET Downsamples a uCT dataset. 
%
%   base_dir is $HOME/Documents/phd/ and set in utils/baseDir()
%
%   Input:
%    - dir_path, path to the directory containing the dataset from base_dir
%    - base_name, name of the dataset.
%    - new_resolution, either an array with the X, Y, and Z resolutions or
%    the factor by which to downsample by.
%    - batch_size, number of image to process in one batch, default 512.
%    - img_extension, extension of the images in the uCT dataset, default
%    bmp.
%   - save_extension, extension of the downsampled images, default png.
%
%   Return:
% Author Alys Clark
% Modified by Emily-Jade Yee
% Modified by Mathias Roesler Nov 2022
narginchk(3, 6);

new_resolution = varargin{1};

if nargin < 6
    save_extension = "png";

else
    save_extension = varargin{4};
end

if nargin < 5
    img_extension = "bmp";

else
    img_extension = varargin{3};
end

if nargin < 4
    batch_size = 512;

else
    batch_size = varargin{2};
end

% Directory where images are located
load_directory = join([baseDir(), dir_path, base_name], '/');

save_directory = join([load_directory, "downsampled"], '/');

%% Load and set parameters
% Load properties of the original images that dont change
toml_map = toml.read(join([load_directory, base_name + ".toml"], '/'));
params = toml.map_to_struct(toml_map);
stack_size = params.end_nb - params.start_nb + 1;

log_file = join([save_directory, params.prefix + "_downsampled.log"], '/');

if sum(size(new_resolution)) == 4
    % The input is the new resolution vector with 3 components
    new_nb_pixel_x = new_resolution(1);
    new_nb_pixel_y = new_resolution(2);
    new_stack_size = new_resolution(3);
    z_factor = new_resolution(3);

elseif isscalar(new_resolution)
    % The input is the downsampling factory
    new_nb_pixel_x = round(params.nb_pixel_x/new_resolution(1));
    new_nb_pixel_y = round(params.nb_pixel_y/new_resolution(1));
    new_stack_size = round(stack_size/new_resolution(1));
    z_factor = new_resolution(1);

else
    error("The size of the input resolution is wrong. Size should be 1 or 3.")
end

new_resolution_x = (params.nb_pixel_x/new_nb_pixel_x)*params.resolution; % um/pixel x-dir
new_resolution_y = (params.nb_pixel_y/new_nb_pixel_y)*params.resolution; % um/pixel y-dir
new_resolution_z = (stack_size/new_stack_size)*params.resolution; % um/pixel z-dir

nb_runs = ceil(stack_size/batch_size); % Number of times to run loop
img_paths = getImagePaths(load_directory, img_extension);
img_save_index = 0;

%% Main resizing loop
for run = 1:nb_runs
    disp("Running batch number " + num2str(run) + "/" + num2str(nb_runs));
    first_image_nb = (run-1) * batch_size + 1;
    last_image_nb = run * batch_size + 1;

    if last_image_nb > params.end_nb
        last_image_nb = params.end_nb;
    end

    batch_stack_size = last_image_nb - first_image_nb;
    new_batch_stack_size = round(batch_stack_size / z_factor);

    % Load all images in current batch
    disp("Loading " + num2str(batch_stack_size) + " images in batch");
    img_stack = loadImageStack(img_paths(first_image_nb:last_image_nb));

    % Resize current batch stack
    new_stack = imresize3(uint8(img_stack), [new_nb_pixel_x, ...
        new_nb_pixel_y, new_batch_stack_size]); 

    % Save current batch stack
    disp("Saving " + num2str(new_batch_stack_size) + " downsampled images");
    saveImageStack(new_stack, save_directory, params.prefix, ...
        img_save_index, save_extension);
    img_save_index = img_save_index + new_batch_stack_size;
end

%% Log required information
file_ID = fopen(log_file, 'w');
fprintf(file_ID, "nb_pixel_x=%3d\n", new_nb_pixel_x);
fprintf(file_ID, "nb_pixel_y=%3d\n", new_nb_pixel_y);
fprintf(file_ID, "stack_size=%3d\n", new_stack_size);
fprintf(file_ID, "pixel_x_res=%03f\n", new_resolution_x);
fprintf(file_ID, "pixel_y_res=%03f\n", new_resolution_y);
fprintf(file_ID, "pixel_z_res=%03f\n", new_resolution_z);
fclose(file_ID);