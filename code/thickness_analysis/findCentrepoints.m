function centrepoints = findCentrepoints(mask, region)
%FINDCENTREPOINTS Finds the centrepoint of a mask given the region.
%
%   If there is no clear separation between the left and right horn, three
%   points are given: left, centre, and right. Otherwise a single point
%   corresponding to the desired region is returned.
%
%   Input:
%    - mask, binary mask.
%    - region, either left or right and used to sort the centrepoints. 
%   Return:
%    - centrepoints, coordinates of the centrepoints,
%    centrepoints(1, 2) or centrepoints(3, 2).
if ~islogical(mask)
    % Convert to logical if not because regionprops only works on logical
    mask = imbinarize(mask);
end

filled_mask = imfill(mask, "holes");
centre_regions = and(not(mask), filled_mask);

props = regionprops(filled_mask, 'Area', 'Centroid');
nb_regions = length(props);
proper_region = 0;

% Go through regions to filter out the abherent regions
for k = 1:nb_regions
    if props(k).Area > 250
        proper_region = proper_region + 1;
        region_idx = k; % Only care about this if there is 1 proper region
    end
end

if proper_region == 1
    % No clear separation between left and right horn, need 3 points
    % Or single horn, need 1 point
    area = props(region_idx).Area;

    % Use area to sort if in the body or if single horn situation
    % 3500 may be dataset depend and a better way of doing might be needed
    if area >= 3500
        % In the body
        centrepoints = zeros(3, 2);
        centrepoints(2, :) = props(region_idx).Centroid; 

    else
        % Single horn and can exit early
        centrepoints = props(region_idx).Centroid;
        return;
    end

else
    centrepoints = zeros(1, 2);
end

props = regionprops(centre_regions, 'Area', 'Centroid');
nb_regions = length(props);

% Create arrays to recuperate region properties
areas = zeros(nb_regions, 1);
centroids = zeros(nb_regions, 2);

for k = 1:nb_regions
    areas(k) = props(k).Area;
    centroids(k, :) = props(k).Centroid;
end

[~, ind] = maxk(areas, 2); % Find two largest areas indices
centroids = centroids(ind, :); % Get the correct centroids

if size(centrepoints, 1) == 3
    % Find the left and right centroids and sort them
    if centroids(1, 1) > centroids(2, 1)
        centrepoints(1, :) = centroids(2, :);
        centrepoints(3, :) = centroids(1, :);

    else
        centrepoints(1, :) = centroids(1, :);
        centrepoints(3, :) = centroids(2, :);
    end

else
    if strcmp(region, "left")
        % Get the left most centre point
        [~, idx] = min(centroids(:, 1));
        centrepoints(1, :) = centroids(idx, :);

    elseif strcmp(region, "right")
        % Get the right most centre point
        [~, idx] = max(centroids(:, 1));
        centrepoints(1, :) = centroids(idx, :);
    end
end

end