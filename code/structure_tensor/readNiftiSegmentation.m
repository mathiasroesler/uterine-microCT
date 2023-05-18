function [img_stack, mask_stack] = readNiftiSegmentation(base_name)
%READNIFTISEGMENTATION Reads the segmentation and images from the nifti
%files. The masks are applied to the images.
%   Input:
%    - base_name, name of the dataset.
%
%   Return:
%    - img_stack, stack of images read from the nifti file. The
%    segmentation masks have been applied.
%    - mask_stack, stack of segmentation masks.
base_dir = join([getenv("HOME"), "Documents/phd/microCT/data"], '/');
src_dir = join([base_dir, base_name, "downsampled/"], '/');
mask_file = base_name + "_segmentation.nii.gz";
img_file = base_name + ".nii.gz";

mask_stack = niftiread(src_dir + mask_file);
img_stack = niftiread(src_dir + img_file);

% Binarize the masks and convert to uint8
mask_stack = uint8(imbinarize(mask_stack)); 

% Permute the first two columns of images and masks
mask_stack = permute(mask_stack, [2 1 3]);
img_stack = permute(img_stack, [2, 1, 3]);

% Apply masks
img_stack = mask_stack .* img_stack;

end