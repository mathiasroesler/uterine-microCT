function centrepoints = findCentrepoints(mask)
%FINDCENTREPOINTS Finds the centrepoints of a mask stack
%
%   Input:
%    - mask, binary mask.
%
%   Return:
%    - centrepoints, centrepoints of the different regions in the mask.
%   The first row is the left region and the second row is the right
%   region.
filled_mask = imfill(mask, "holes");
centre_region = and(not(mask), filled_mask);
props = regionprops(centre_region, 'Area', 'Centroid');

% Create arrays to recuperate region properties
areas = zeros(3, 1);
centroids = zeros(3, 2);

for k = 1:length(props)
    areas(k) = props(k).Area;
    centroids(k, :) = props(k).Centroid;
end

[~, ind] = maxk(areas, 2); % Find two largest areas indices
centrepoints = centroids(ind, :); % Get the correct centroids

if centrepoints(1) > centrepoints(2)
    % Swap the centrepoints to ensure the left one is first
    centrepoints([1, 2], :) = centrepoints([2, 1], :);
end