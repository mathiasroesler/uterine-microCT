function downsampleMicroCTDataset(dir_name, primary_direction_size, ...
    new_stack_size, batch_size, img_extension, save_extension)
%DOWNSAMPLEMICROCTDATASET Downsamples a uCT dataset. 
%
%   Input:
%    - dir_name, name of the uCT dataset directory.
%    - primary_direction_size, new size of the longest direction of the
%    image. The aspect ration is maintained so only one direction is
%    specified.
%    - new_stack_size, new number of slices in the image stack.
%    - batch_size, number of image to process in one batch, default 512.
%    - img_extension, extension of the images in the uCT dataset, default
%    bmp.
%   - save_extension, extension of the downsampled images, default png.
%
%   Return:
% Author Alys Clark
% Modified by Emily-Jade Yee
% Modified by Mathias Roesler Nov 2022
if nargin < 6
    save_extension = "png";
end

if nargin < 5
    img_extension = "bmp";
end

if nargin < 4
    batch_size = 512;
end

load_directory = join([getenv("HOME"), ...
    "Documents/phd", dir_name], ...
    '/'); % Directory where images are located
save_directory = join([load_directory, "downsampled"], '/');

%% Load and set parameters
% Load properties of the original images that dont change
[nb_pixel_x, nb_pixel_y, resolution, start_nb, end_nb, img_prefix] ...
    = loadParams(load_directory);
stack_size = end_nb - start_nb;

log_file = join([save_directory, img_prefix + "_downsampled.log"], '/');

if nb_pixel_x >= nb_pixel_y
    aspect_ratio = nb_pixel_y / nb_pixel_x;
    new_nb_pixel_x = primary_direction_size;
    new_nb_pixel_y = round(aspect_ratio * primary_direction_size);

else
    aspect_ratio = nb_pixel_x / nb_pixel_y;
    new_nb_pixel_y = primary_direction_size;
    new_nb_pixel_x = round(aspect_ratio * primary_direction_size);
end

new_resolution_x = (nb_pixel_x/new_nb_pixel_x)*resolution; % um/pixel x-dir
new_resolution_y = (nb_pixel_y/new_nb_pixel_y)*resolution; % um/pixel y-dir
new_resolution_z = (stack_size/new_stack_size)*resolution; % um/pixel z-dir

nb_runs = ceil(stack_size/batch_size); % Number of times to run loop

% Get the start and end numbers of the images in each batch
first_image_nbs = start_nb + ((batch_size+1) * (0:nb_runs-1));
last_image_nbs = batch_size + first_image_nbs;
last_image_nbs(end) = end_nb; % Reset the last number for the last batch
nb_img_per_batch = last_image_nbs - first_image_nbs;
tmp_stack_sizes = round(((nb_img_per_batch) * new_stack_size) ...
    / stack_size); % Stack size for each batch for resize

if sum(tmp_stack_sizes) ~= new_stack_size
    tmp_stack_sizes(end) = tmp_stack_sizes(end) + new_stack_size - sum( ...
        tmp_stack_sizes);
end

%% Main resizing loop
for run = 1:nb_runs
    disp("Running batch number " + num2str(run) + "/" + num2str(nb_runs));
    batch_stack_size = tmp_stack_sizes(run);
    img_stack = zeros(nb_pixel_x, nb_pixel_y, batch_stack_size);
    save_img_cnt = (1:batch_stack_size) + tmp_stack_sizes(1) * (run-1);
    img_nb_list = first_image_nbs(run):last_image_nbs(run);

    % Load and stack all images in a batch
    disp("Loading " + num2str(nb_img_per_batch(run)) + " images in batch");
    parfor j = 1:length(img_nb_list)
        img_nb = join([sprintf("%08.f", img_nb_list(j)), ...
            img_extension], '.');
        img_name = join([img_prefix, img_nb], '_');

        img_stack(:, :, j) = imread(join([load_directory, ...
            img_name], '/'));
    end

    new_stack = imresize3(uint8(img_stack), [new_nb_pixel_x, ...
        new_nb_pixel_y, batch_stack_size]); % Resize current batch stack
    
    disp("Saving " + num2str(batch_stack_size) + " downsampled images");
    parfor k = 1:batch_stack_size
        img_nb = join([sprintf("%03.f", save_img_cnt(k)), ...
            save_extension], '.');
        img_name = join([img_prefix, img_nb], '_');

        imwrite(new_stack(: ,:, k), join([save_directory, img_name], ...
            '/'), save_extension);
    end
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