function [img_stack, mask_stack] = readNiftiSegmentation(dir_path, ...
    base_name, downsampled)
%READNIFTISEGMENTATION Reads the segmentation and images from the nifti
%files. The masks are applied to the images.
%   
%   base_dir is $HOME/Documents/phd/
%
%   Input:
%    - dir_path, path to the directory containing the dataset from base_dir
%    - base_name, name of the dataset.
%    - downsampled, true if the dataset has been downsampled, default value
%    is true.
%
%   Return:
%    - img_stack, uint8, stack of images read from the nifti file. The
%    segmentation masks have been applied.
%    - mask_stack, logical, stack of segmentation masks.
if nargin < 3
    downsampled = true;
end

load_directory = join([getenv("HOME"), "Documents/phd", dir_path, base_name], '/');

if downsampled
    load_directory = join([load_directory, "downsampled"], '/');
end

mask_file = base_name + "_segmentation.nii.gz";
img_file = base_name + ".nii.gz";

mask_stack = niftiread(join([load_directory, mask_file], '/'));
img_stack = niftiread(join([load_directory, img_file], '/'));

% Binarize the masks
mask_stack = imbinarize(mask_stack); 

% Permute the first two columns of images and masks
mask_stack = permute(mask_stack, [2 1 3]);
img_stack = permute(img_stack, [2, 1, 3]);

% Apply masks
img_stack = uint8(mask_stack) .* img_stack;

end