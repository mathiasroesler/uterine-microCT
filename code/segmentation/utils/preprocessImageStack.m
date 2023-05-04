function processed_img_stack = preprocessImageStack(img_stack, fct_nb, ...
    varargin)
%PREPROCESSIMAGE Preprocesses the inputed image based on the selected
%method.
%   Inputs:
%    - img, image to preprocess.
%    - fct_nb, number associated with the preprocessing function:
%     1 -> imadjust
%     2 -> histeq
%     3 -> imsharpen
%     4 -> imflatfield
%    - varargin, extra arguments for the imflatfield function.
narginchk(2, 4);
sigma = 0;
mask_stack = zeros(size(img_stack));

if nargin == 3
    img_strel = strel('disk', varargin{1});

elseif nargin == 4
    sigma = varargin{1};
    mask_stack = varargin{2};
end

processed_img_stack = uint8(zeros(size(img_stack)));

parfor k = 1:size(processed_img_stack, 3)
    switch fct_nb
        case 1
            processed_img_stack(:, :, k) = imadjust(img_stack(:, :, k));
            processed_img_stack(:, :, k) = imclose( ...
                processed_img_stack(:, :, k), img_strel);
        case 2
            processed_img_stack(:, :, k) = histeq(img_stack(:, :, k));
            processed_img_stack(:, :, k) = imclose( ...
                processed_img_stack(:, :, k), img_strel);
        case 3
            processed_img_stack(:, :, k) = imsharpen(img_stack(:, :, k));
            processed_img_stack(:, :, k) = imclose( ...
                processed_img_stack(:, :, k), img_strel);    
        case 4
            processed_img_stack(:, :, k) = imflatfield(img_stack(:, :, k), ...
                sigma, mask_stack(:, :, k));
    end

end