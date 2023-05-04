function saveImageStack(img_stack, save_dir, base_name, start_nb, ...
    extension)
%SAVEIMAGESTACK Writes the images in the stack to the save directory.
%   
%   Inputs:
%    - img_stack, images to save.
%    - save_dir, path to the directory in which to save the images.
%    - base_name, base name for the images.
%    - start_nb, number at which to start saving the stack, default 1.
%    - extension, image extension, default pmg.
if nargin < 5
    extension = "png";
end

if nargin < 4
    start_nb = 1;
end

parfor k = 1:size(img_stack, 3)
    img_name = sprintf(base_name + "_%.3d." + extension, (k-1) + start_nb);
    imwrite(img_stack(:, :, k), join([save_dir, img_name], '/'));
end