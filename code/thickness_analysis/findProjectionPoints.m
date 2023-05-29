function projection_points = findProjectionPoints(img, centre_point, nb_points)
    % Find the projection points from the centre point onto the muscle
    % layers in four directions (top, bottom, left, right)

    % Check if the number of points is even
    if mod(nb_points, 2) ~= 0
        error('Error: the number of points needs to be even.');
    end

    % Ensure that coordinates are integers
    x_centre = round(centre_point(2));
    y_centre = round(centre_point(1));

    angles = linspace(1, 1 + (nb_points/2), nb_points/2) * (pi / (nb_points / 2));
    projection_points = zeros(nb_points*2, 2);

    for i = 1:numel(angles)
        theta = angles(i);

        if theta ~= pi
            y_points = (1:img.shape(1))';
            x_points = (((y_centre - y_points) * cos(theta)) / sin(theta)) + x_centre;
            x_points = int32(x_points); % Convert to int

            % Find the points that are in the image
            intersection = intersect(find(x_points >= 0), find(x_points < size(img, 2)));
            x_points = x_points(intersection);
            y_points = y_points(intersection);

        else
            x_points = (1:img.shape(2))';
            y_points = ones(size(x_points)) * y_centre;
        end

        line_indices = skd.line(y_points(1), x_points(1), y_points(end), x_points(end));
        line_y = line_indices{1};
        line_x = line_indices{2};
        line = img(sub2ind(size(img), line_y, line_x));
        coords = find(line(1:end-1) ~= line(2:end)) + 1;

        for j = 1:2:numel(coords)
            coords(j) = coords(j) + 1;
        end

        projection_points(((i-1)*4)+1:i*4, :) = createProjectionPointCoords(...
            line_x(coords), line_y(coords), centre_point);
    end
end
