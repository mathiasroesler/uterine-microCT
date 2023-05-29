function [muscle_thickness_array, slice_thickness_array] = estimateMuscleThickness(img_stack, centreline, nb_points, slice_nbs)
    % Estimates the muscle thickness of each slice

    nb_imgs = size(img_stack, 1);
    muscle_thickness_array = zeros(nb_imgs, 1);
    slice_thickness_array = cell(1, numel(slice_nbs));

    for i = 1:nb_imgs
        img = squeeze(img_stack(i, :, :));
        projection_points = findProjectionPoints(img, centreline(i, :), nb_points);

        diff = diff(projection_points);
        norm = vecnorm(diff, 2, 2);
        thickness = norm(1:2:end);
        muscle_thickness_array(i) = mean(thickness);

        if ismember(i, slice_nbs)
            % Order thickness to go from 0 to 2pi
            ordered_thickness = [thickness(1:2:end); thickness(2:2:end)];
            ordered_thickness = ordered_thickness(:);

            % Roll array to line up 0 with anti-mesometrial border
            [~, max_idx] = max(ordered_thickness);
            ordered_thickness = circshift(ordered_thickness, nb_points - max_idx);
            slice_thickness_array{slice_nbs == i} = ordered_thickness;
        end
    end

    slice_thickness_array = cat(2, slice_thickness_array{:});
end
