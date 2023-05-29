function projection_points = createProjectionPointCoords(x_coords, y_coords, centre_point)
%CREATEPROJECTIONPOINTCOORDS  Creates the (x, y) pairs of coordinates for 
% the projection points
% 
% 	Input:
% 	 - x_coords, list of coordinates of the projection points on the x axis.
% 	 - y_coords, list of coordinates of the projection points on the y axis.
% 	 - centre_point, coordinates of the centre point.
% 	
% 	Return:
% 	 - projection_points, list of coordinates of the projection points.

% Check if x_coords and y_coords have the same size
if numel(x_coords) ~= numel(y_coords)
    error('Error: x_coords and y_coords should have the same size.');
end

point_list = [y_coords(:), x_coords(:)];

if numel(x_coords) == 4
    projection_points = point_list;
    return;
end

diff = point_list - centre_point;
points_below = find(diff(:, 2) < 0); % Before x centre
points_above = find(diff(:, 2) > 0); % After x centre

if isempty(points_below)
    points_below = find(diff(:, 1) > 0);
    points_above = find(diff(:, 1) < 0);
end

[~, lower_idx] = min(vecnorm(diff(points_below, :), 2, 2));
[~, upper_idx] = min(vecnorm(diff(points_above, :), 2, 2));

if diff(upper_idx, 1) < 0
    upper_points = point_list(upper_idx-1:upper_idx, :);
    lower_points = point_list(lower_idx:lower_idx+1, :);
    projection_points = [upper_points; lower_points];
elseif diff(upper_idx, 1) >= 0
    upper_points = point_list(upper_idx:upper_idx+1, :);
    lower_points = point_list(lower_idx-1:lower_idx, :);
    projection_points = [lower_points; upper_points];
end
end
