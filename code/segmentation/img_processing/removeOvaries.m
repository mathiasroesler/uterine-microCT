function removeOvaries(dir_name, left_start_nb, right_start_nb, ...
    params_file, extension, start_nb)
%REMOVEOVARIES Removes the ovaries from the segmentation masks of a uCT
%dataset.
%
%   Input:
%    - dir_name, name of the uCT dataset directory.
%    - left_start_nb, image in the stack after which to remove the left
%    ovary.
%    - right_start_nb, image in the stack after which to remove the right
%    ovary.
%    - params_file, name of the parameters file to read from, default
%    preprocess.params.
%    - extension, extension of the segmentation masks, default png.
%    - start_nb, number at which to start saving images, default 1.
%   Return:
if nargin < 6
    start_nb = 1;
end

if nargin < 5
    extension = "png";
end

if nargin < 4
    params_file = "preprocess.params";
else
    [~, ~, extension] = fileparts(params_file);

    % Add the file name if not provided
    if extension ~= ".params"
        params_file = params_file + ".params";
    end
end

dir_path = join([getenv("HOME"), "Documents/phd", dir_name], '/');
img_paths = getImagePaths(dir_path, extension);

disp('Loading image stack')
img_stack = loadImageStack(img_paths);

middle_pixel = floor(max(size(img_stack)) / 2);
[~, largest_dim] = max(size(img_stack));
nb_img = size(img_stack, 3);

disp('Removing left ovary')
for k = left_start_nb:nb_img
    % Remove ovary on the left side of the image
    switch largest_dim
        case 1
            img_stack(1:middle_pixel, :, k) = 0;

        case 2
            img_stack(:, 1:middle_pixel, k) = 0;
    end
end

disp('Removing right ovary')
for k = right_start_nb:nb_img
    % Remove ovary on the right side of the image
    switch largest_dim
        case 1
            img_stack(middle_pixel:end, :, k) = 0;

        case 2
            img_stack(:, middle_pixel:end, k) = 0;
    end
end

disp('Saving image stack')
img_prefix = loadParams(dir_path + '/' + params_file);
saveImageStack(img_stack, dir_path, img_prefix, start_nb, extension);

end