function centrepoints = findCentrepoints(mask, horn)
%FINDCENTREPOINTS Finds the centrepoints of a mask stack and places the one
%corresponding to the horn first
%
%   Input:
%    - mask, binary mask.
%    - horn, either left or right. Finds
%   Return:
%    - centrepoints, centrepoints of the different regions in the mask.
%   The first row is the centre point that corresponds to the inputed horn.
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

if centrepoints(1) > centrepoints(2) && strcmp(horn, "left") && centrepoints(2) > 0
    % Swap the centrepoints to ensure the left one is first
    % Make sure the second centre point exists first
    centrepoints([1, 2], :) = centrepoints([2, 1], :);

elseif centrepoints(2) > centrepoints(1) && strcmp(horn, "right")
    % Swap the centrepoints to ensure the right one is first
    centrepoints([2, 1], :) = centrepoints([1, 2], :);
end

end