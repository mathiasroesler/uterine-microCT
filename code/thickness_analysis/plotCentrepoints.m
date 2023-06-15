dir_name = "AWA015_PTA_1_Rec_Trans/downsampled/muscle_segmentation/right";
extension = "png";
base_dir = join([getenv("HOME"), "Documents/phd/microCT/data", dir_name], '/');

mask_paths = getImagePaths(base_dir, extension);
mask_stack = loadImageStack(mask_paths);

for k = 1:size(mask_stack, 3)
    mask = mask_stack(:, :, k);
    centrepoints = findCentrepoints(mask, "right");

    fig = figure;
    set(gcf,'visible','off')
    imshow(imbinarize(mask)); hold on
    for j = 1:size(centrepoints, 1)
        plot(centrepoints(j, 1), centrepoints(j, 2), '.r')
    end
    print(fig, '-dpng', base_dir + "/test/" + k)    
    close all;
    clear fig;
end