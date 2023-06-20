function angle = ComputeFibreAngle(fibre, centrepoints, cur_X, cur_Z)
%COMPUTEFIBREANGLE Computes the angle in degrees of the fibre relative to
%the plane given by the vector between two centre points. 
%
%   Inputs:
%    - fibre, 3D array containing the XYZ coordinates of the fibre.
%    - centrepoints, 6xN, list of centrepoints.
%    - cur_X, x coordinate of the current fibre in the general coordinate
%    system.
%    - cur_Z, z coordinate of the current fibre in the general coordinate
%    system.
%
%   Return:
%    - angle, angle between the fibre and the current plane in degrees.
if isempty(centrepoints)
    z_vector = [0; 0; 1];

else
    cur_Z = round(cur_Z);
    cur_X = round(cur_X);
    
    if all(centrepoints(3:4, cur_Z))
        % If a middle point is found
        % Define the z vector as the centre vector between current slice
        % and the current + 5th slice.
        cur_idx = 3;
        cur_centrepoint = centrepoints(cur_idx:cur_idx+1, cur_Z);

    elseif ~all(centrepoints(5:6, cur_Z))
        % If there is only one centre point found
                cur_idx = 5;
        cur_centrepoint = centrepoints(cur_idx:cur_idx+1, cur_Z);

    else
        % If there are left and right centre points find the neares to
        % cur_X
        [~, cur_idx] = min(abs(centrepoints(1:2:end, cur_Z) - cur_X));
        cur_centrepoint = centrepoints(cur_idx:cur_idx+1, cur_Z);

    end

    if cur_Z + 5 > size(centrepoints, 2)
        % If there are not enough slices use the previous slices
        % instead
        next_centrepoint = centrepoints(cur_idx:cur_idx+1, cur_Z - 5);
        z_vector = cur_centrepoint - next_centrepoint;
    else
        next_centrepoint = centrepoints(cur_idx:cur_idx+1, cur_Z + 5);
        z_vector = next_centrepoint - cur_centrepoint;
    end
    z_vector = [z_vector; 5]; % Append the z component
end

angle = rad2deg(acos(fibre * z_vector));
end