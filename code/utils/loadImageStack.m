function img_stack = loadImageStack(img_paths, mask_paths)
%LOADIMAGESTACK Loads the images contains in the img_paths cell array.
%
% If mask_paths is specified, the masks are applied to the images.
%   
%   Input:
%    - img_paths, cell array containing the path to every image,
%    img_paths(N x 1).
%    - mask_paths, cell array containing the path to every mask,
%    mask_paths(N x 1), default value [].
%
%   Return:
%    - img_stack, matrix containing the images read from the paths. If
%    mask_paths is specified, the masks are applied before return.
if nargin < 2
    mask_paths = [];
end

nb_imgs = length(img_paths);
if nb_imgs ~= length(mask_paths) && ~isempty(mask_paths)
    error("The number of images is different to the number of masks");
end

img = imread(img_paths{1}); % Read first image for size
img_stack = zeros(size(img, 1), size(img, 2), nb_imgs, 'uint8'); 

if ~isempty(mask_paths)
    for k = 1:nb_imgs
        img = imread(img_paths{k}); 
        mask = imread(mask_paths{k});

        if isnumeric(mask)
            % Binarize the image if it is not of logical type
            mask = imbinarize(mask);
        end

        mask = uint8(mask);
        img_stack(:, :, k) = img .* mask;
    end

else
    for k = 1:nb_imgs
        img_stack(:, :, k) = imread(img_paths{k});
    end

end