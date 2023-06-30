function centrepoints = findCentrepoints(mask, region)
%FINDCENTREPOINTS Finds the centrepoint of a mask given the region.
%
%   If there is no clear separation between the left and right horn, three
%   points are given: left, centre, and right.
%
%   Input:
%    - mask, binary mask.
%    - region, place the centrepoint in the corresponding region if only
%    one centrepoint is found.
%
%   Return:
%    - centrepoints, coordinates of the centrepoints, the first pair of
%    is the left point, the second pair the middle point, the third pair
%    the right point, centrepoints(3, 2).
if ~islogical(mask)
    % Convert to logical if not because regionprops only works on logical
    mask = imbinarize(mask);
end

proper_size = 300; % Define how big a region is to be considered
centrepoints = zeros(3, 2); % Placeholder for the centrepoints

filled_mask = imfill(mask, "holes");
centre_regions = and(not(mask), filled_mask);

filled_props = regionprops(filled_mask, 'Area', 'Centroid');
proper_region = sum([filled_props.Area] > proper_size);

centre_props = regionprops(centre_regions, 'Area', 'Centroid');
nb_holes = sum([centre_props.Area] > proper_size);

if proper_region == 1
    % No clear separation between left and right horn, need 3 points
    % Or single horn, need 1 point
    % Use the number of holes to filter out the region. If there are two
    % then in the body, if there is 1 then single horn
    if nb_holes > 1
        % In the body
        centrepoints(2, :) = filled_props( ...
            [filled_props.Area] > proper_size).Centroid; 

    else
        % Single horn and can exit early
        centrepoint = filled_props( ...
            [filled_props.Area] > proper_size).Centroid;

        if strcmp(region, "left")
            % Left horn
            centrepoints(1, :) = centrepoint;

        else
            % Right horn
            centrepoints(3, :) = centrepoint;
        end

        return;
    end
end

% Create arrays to recuperate region properties if not exited early
areas = zeros(length(centre_props), 1);
centroids = zeros(length(centre_props), 2);

for k = 1:length(centre_props)
    areas(k) = centre_props(k).Area;
    centroids(k, :) = centre_props(k).Centroid;
end

[~, ind] = maxk(areas, 2); % Find two largest areas indices
centroids = centroids(ind, :); % Get the correct centroids

% Refine the centres by using the skeleton
skeleton = bwskel(centre_regions);
[idx_y, idx_x] = find(skeleton == 1);

for k = 1:size(centroids, 1)
    differences = [idx_x, idx_y] - centroids(k, :);
    [~, min_idx] = min(vecnorm(differences'));
    centroids(k, :) = [idx_x(min_idx), idx_y(min_idx)];
end

% Find the left and right centroids and sort them
if centroids(1, 1) > centroids(2, 1)
    centrepoints(1, :) = centroids(2, :);
    centrepoints(3, :) = centroids(1, :);

else
    centrepoints(1, :) = centroids(1, :);
    centrepoints(3, :) = centroids(2, :);
end

end