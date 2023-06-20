function ratio = estimateMuscleRatio(dir_name, downsampled, extension)
%ESTIMATEMUSCLERATIO Estimates the quantity of myometrium versus
%endometrium based on the tissue and muscle segmentation masks
%
%   Input:
%    - dir_name, name of the uCT dataset directory.
%    - downsampled, 1 if the uCT dataset was downsampled first and 0
%       otherwise, default 1.
%    - extension, extension of the segmentation masks, default png.
%
%   Return:
%    - ratio, muscle to tissue ratio for each slice.
if nargin < 2
    downsampled = 1;
end
if nargin < 3
    extension = "png";
end

if downsampled
    muscle_dir = "downsampled/muscle_segmentation";
    tissue_dir = "downsampled/tissue_segmentation";
else
    muscle_dir = "muscle_segmentation";
    tissue_dir = "tissue_segmentation";
end

muscle_dir_path = join([getenv("HOME"), "Documents/phd", dir_name, ...
    muscle_dir], '/');
tissue_dir_path = join([getenv("HOME"), "Documents/phd", dir_name, ...
    tissue_dir], '/');
muscle_mask_paths = getImagePaths(muscle_dir_path, extension);
tissue_mask_paths = getImagePaths(tissue_dir_path, extension);


% Load masks and convert to logicals
muscle_mask_stack = logical(loadImageStack(muscle_mask_paths));
tissue_mask_stack = logical(loadImageStack(tissue_mask_paths));

% Size checking
if ~ all(size(muscle_mask_stack) == size(tissue_mask_stack))
    error("The size of the tissue and muscle segmentation mask do not agree");
end

muscle_stack_nb_pixels = squeeze(sum(muscle_mask_stack, [1, 2]));
tissue_stack_nb_pixels = squeeze(sum(tissue_mask_stack, [1, 2]));
ratio = muscle_stack_nb_pixels ./ tissue_stack_nb_pixels;

end