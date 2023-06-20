function segmentMicroCTDataset(dir_name, segmentation_type, split_nb, ...
    extension, start_nb)
%SEGMENTMICROCTDATASET Segments a uCT dataset based on the segmentation 
%type. 
%
%   Input:
%    - dir_name, name of the uCT dataset directory.
%    - segmentation_type, type of segmentation to use:
%      1 is tissue segmentation
%      2 is muscle layer segmentation
%      3 is fat segmentation
%      4 is shape segmentation
%    - split_nb, slice at which to split into left and right horn.
%    - extension, extension of the segmentation masks, default png.
%    - start_nb, number at which to start saving images, default 1.
%   Return:
if nargin < 5
    start_nb = 1;
end
if nargin < 4
    extension = "png";
end

dir_path = join([getenv("HOME"), "Documents/phd", dir_name], '/');
img_paths = getImagePaths(dir_path, extension);
[img_prefix, preprocess, morph_size, sigma] = loadParams(dir_path + ...
    "/preprocess.params");

disp('Loading image stack')
img_stack = loadImageStack(img_paths);

disp('Spliting image stack')
stack_cell = splitImageStack(img_stack, split_nb);
mask_cell = cell(size(stack_cell)); % Empty cell array for the masks

for k = 1:length(stack_cell)
    img_stack = stack_cell{k};
    disp('Processing stack ' + string(k))
    if preprocess ~= 0
        disp('    Preprocessing image stack')
        if preprocess == 4
            mask_paths = getImagePaths(dir_path + "/tissue_segmentation", ...
                extenstion);
            mask_stack = loadImageStack(mask_paths);
            img_stack = preprocessImageStack(img_stack, preprocess, sigma, ...
                mask_stack);
            clear mask_paths mask_stack;
        else
            img_stack = preprocessImageStack(img_stack, preprocess, ...
                morph_size);
        end
    end

    disp('    Segmenting image stack')
    switch segmentation_type
        case 1
            % Tissue segmentation
            save_dir= join([dir_path, "tissue_segmentation"], '/');
            [~, threshold, morph_size, neighborhood_size] = loadParams( ...
                save_dir + "/segmentation.params");
            mask_stack = tissueSegmentation(img_stack, threshold, morph_size, ...
                neighborhood_size);

        case 2
            % Muscle segmentation
            save_dir = join([dir_path, "muscle_segmentation"], '/');
            [~, morph_size, neighborhood_size] = loadParams(save_dir + ...
                "/segmentation.params");
            mask_stack = muscleSegmentation(img_stack, morph_size, ...
                neighborhood_size);

        case 3
            % Fat segmentation
            save_dir= join([dir_path, "fat_segmentation"], '/');
            [~, morph_size] = loadParams(save_dir + "/segmentation.params");
            mask_stack = fatSegmentation(img_stack, morph_size);

        case 4
            % Shape segmentation
            save_dir= join([dir_path, "shape_segmentation"], '/');
            [~, threshold, morph_size, neighborhood_size] = loadParams( ...
                save_dir + "/segmentation.params");
            mask_stack = shapeSegmentation(img_stack, threshold, morph_size, ...
                neighborhood_size);
    end

    mask_cell{k} = mask_stack;
end

disp('Fusing mask stack')
mask_stack = fuseImageStacks(mask_cell);

disp('Saving mask stack')
saveImageStack(mask_stack, save_dir, img_prefix, start_nb, extension);


end