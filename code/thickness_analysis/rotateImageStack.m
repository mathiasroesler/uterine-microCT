function rotated_stack = rotateImageStack(img_stack, horn, nb_used_slices)
%ROTATEIMAGESTACK Rotates each image in the stack based on the direction of
%the centreline vector
%
%   Input:
%    - img_stack, stack of N MxP images to rotate, img_stack(MxPxN).
%    - horn, selects which horn to rotate, either left or right.
%    - nb_used_slices, number of slices to take to estimate the centre
%    vector.
%
%   Return:
%    - rotated_stack, stack of N-nb_used_slices MxP rotated images, 
%   rotated_stack(MxPx(N-nb_used_slices)).
rotated_stack = zeros(size(img_stack));

if matches(horn, 'left')
    direction = 1;
    
elseif matches(horn, 'right')
    direction = 2;
else
    error("Error: incorrect value for horn. Should be either left or right")
end

if ~isnumeric(img_stack)
    % Conver to double if not already numeric
    img_stack = double(img_stack);
end

for k = 1:size(img_stack, 3)-nb_used_slices
    cur_mask = img_stack(:, :, k); % Current mask to rotate
    next_mask = img_stack(:, :, k+nb_used_slices); % Next mask for finding rotation axis

    cur_centrepoints = findCentrepoints(cur_mask);
    next_centrepoints = findCentrepoints(next_mask);

    % Get the normalised centre vector in 3D
    centre_vector = [next_centrepoints(direction, :) - cur_centrepoints(direction, :), nb_used_slices];
    centre_vector = centre_vector ./ norm(centre_vector);

    % Create the transformation matrix
    origin = [cur_centrepoints(direction, :), 0];
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

    rotated_stack(:, :, k) = rotated_mask;
end