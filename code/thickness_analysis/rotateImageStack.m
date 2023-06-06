function rotated_stack = rotateImageStack(img_stack, region, nb_used_slices)
%ROTATEIMAGESTACK Rotates each image in the stack based on the direction of
%the centreline vector
%
%   Input:
%    - img_stack, stack of N MxP images to rotate, img_stack(MxPxN).
%    - region, selects which region to rotate, either left, right or body.
%    - nb_used_slices, number of slices to take to estimate the centre
%    vector.
%
%   Return:
%    - rotated_stack, stack of N MxP rotated images, 
%   rotated_stack(MxPxN).
nb_slices = size(img_stack, 3);
rotated_stack = zeros(size(img_stack));
idx_to_remove = [];

if ~matches(region, 'left') && ~matches(region, 'right') && ~matches(region, 'body')
    error("Error: incorrect value for horn. Should be either left, right or body")
end

if ~isnumeric(img_stack)
    % Conver to double if not already numeric
    img_stack = double(img_stack);
end

for k = 1:nb_slices
    cur_mask = img_stack(:, :, k); % Current mask to rotate

    % Next mask for finding rotation axis
    if k < nb_slices-nb_used_slices
        next_mask = img_stack(:, :, k+nb_used_slices); 
    else
        % Use the previous slices to get rotation vector
        next_mask = img_stack(:, :, k-nb_used_slices); 
    end
    % Find centre points
    cur_centrepoints = findCentrepoints(cur_mask, region);
    next_centrepoints = findCentrepoints(next_mask, region);

    % Get the normalised centre vector in 3D
    centre_vector = [next_centrepoints(1, :) - cur_centrepoints(1, :), nb_used_slices];
    centre_vector = centre_vector ./ norm(centre_vector);

    % Create the transformation matrix
    origin = [cur_centrepoints(1, :), 0];
    T = findRotationMatrix(centre_vector, [0, 0, 1], origin);

    % Pad the current mask to make a 3D object for rotation
    mask_3D = padarray(cur_mask, [0, 0, 1], 0, 'post');
    
    % Rotate image and collapse it onto the XY plane
    rotated_mask = sum(imwarp(mask_3D, affine3d(T)), 3);

    % Resize the rotated image to have dimensions of original
    dim_diff = size(cur_mask) - size(rotated_mask);
    
    if dim_diff(1) > 0
        rotated_mask = padarray(rotated_mask, [dim_diff(1), 0], 'post');

    else
        rotated_mask = rotated_mask(1:end+dim_diff(1), :);
    end

    if dim_diff(2) > 0
        rotated_mask = padarray(rotated_mask, [0, dim_diff(2)], 'post');

    else
        rotated_mask = rotated_mask(:, 1:end+dim_diff(2));
    end   

    % Check the ratio of white pixel to remove ill-rotated slices
    nb_w_pixels = sum(sum(rotated_mask > 0));

    if nb_w_pixels / numel(rotated_mask) > 1e-3
        % Clean up mask and add rotated mask to stack
        rotated_stack(:, :, k) = imclose(rotated_mask, strel("disk", 1, 4));

    else
        disp("Slice " + k + ": bad rotation");
        idx_to_remove = [idx_to_remove, k];
    end
end

rotated_stack(:, :, idx_to_remove) = [];