
function segtrack_multiplechannels(dir_name, bin, maxfiles, varargin)
    
    % GLOBAL PARAMETERS
    global processnb imx imy  dt_max it iteration at mag thres_offset
    
    % blur window for background estimation
    BLUR_WIN = 101;
    
    % number of processes to run in parallel
    processnb = 2; 
    % image sizes
    imx = 1024/bin; imy = 1024/bin;
    % threshold on intensity
    it = 0.3;
    % distance threshold
    dt_max = 60/bin;
    % area threshold
    at = 0.4;
    % magnification, relative to 20x
    mag = 1;
    
    %% CHECK THIS PARAMETER TRY: 0.005 - 0.02
    thres_offset = 0.05;
    %%
    
    % need three inputs
    if nargin < 4
        error('Use segtrack(directory_name, binning, maxfiles, background_color, channels_to_segment) with full path to the directory of images');
    end
    
    if exist(dir_name, 'dir')
        cd(dir_name);
    else
        error('Invalid directory. Please try again');
    end
    
    % output masks and lineage information written to new folders 
    if ~exist('segmented', 'dir')
        mkdir('segmented');
        mkdir('cellproperties');
    end
    
    if ~exist('matrices', 'dir')
        mkdir('matrices');
    end
    
    numchannels = nargin - 3;
    inputfilenames = cell(1, numchannels);
    for j = 1:numchannels
        channeltosegment = varargin{j};
        
        imagefilenames = import_images(dir_name, channeltosegment);
        
        imagefilenames = imagefilenames(1:maxfiles);
        inputfilenames{1, j} = transpose(imagefilenames);
    end

    % initialize loop variables
    iteration = 1; tot_iterations = floor(maxfiles - 1)/processnb;

    
    im0 = zeros(imx, imy, 'uint16');
    for k = 1:numchannels
        % initialize by segmenting the first image
        % make a 'roster' of cell objects and their properties and
        % write to file
        imagefilenames = inputfilenames{1, k};
        im0 = im0 + imread(imagefilenames{1});
    end

    big_smooth_filt = fspecial('gaussian', floor(BLUR_WIN/bin), 0.3*(BLUR_WIN/(2*bin) - 1) + 0.8);
    min_filt = ordfilt2(double(im0), 10, ones(50), 'symmetric');
    min_filt = imfilter(min_filt, big_smooth_filt, 'replicate', 'conv');
    back_norm_im = double(im0)./min_filt;


    segim0 = laplaceimagesegment_weak(back_norm_im);
    segim0(1:2, 1:2) = 1;
    [segim0, nbcell0] = bwlabel(segim0, 4);
    cell_roster0 = regionprops(segim0, im0, 'Area', 'BoundingBox', 'Centroid', 'Image', 'Orientation', 'Eccentricity', 'WeightedCentroid', 'MajorAxisLength');
    back_vals0 = zeros(length(cell_roster0), 1); centroids = {cell_roster0.Centroid};
    for k = 1:length(cell_roster0)
        back_vals0(k) = find_background(centroids{k}, im0, segim0, bin);
    end
       
    imagefilename = inputfilenames{1}{1};
    outfilename = ['segmented' filesep imagefilename(1:end - 4) '_regionprops.mat'];
    save(outfilename, 'cell_roster0');
    
    % initialize lineage matrix
    lineages = zeros(nbcell0, length(imagefilenames));
    lineages(1:end, 1) = 1:nbcell0;
        

    carryover_im = im0;
    carryover_segim = segim0;
       
    % the main loop
    
    while iteration <= tot_iterations  
        
        i = processnb*(iteration - 1) + 1;

        % read in 'processnb' images and segment them in parallel
        ims = zeros(imx, imy, processnb + 1, 'uint16');
        segims = ims;
        
        ims(:, :, 1) = carryover_im; segims(:, :, 1) = carryover_segim;
        % re-evaluate background illumination every 6 time-points
%         if mod(iteration, 3) == 1 && iteration > 1
%             im0 = carryover_im;
%             big_smooth_filt = fspecial('gaussian', floor(BLUR_WIN/bin), 0.3*(BLUR_WIN/(2*bin) - 1) + 0.8);
%             min_filt = ordfilt2(double(im0), 10, ones(50), 'symmetric');
%             min_filt = imfilter(min_filt, big_smooth_filt, 'replicate', 'conv');
%         end
                  
        % calculates relative shift between two images
        shifts = zeros(processnb, 2);
        
        for k = 1:processnb
            for j = 1:numchannels
                imagefilenames = inputfilenames{1, j};
                ims(:, :, k + 1) = ims(:, :, k + 1) + imread(imagefilenames{i + k});
            end 

            segim1 = laplaceimagesegment_weak(double(ims(:, :, k + 1))./min_filt);
            segim1(1:2, 1:2) = 1;
            [segim1, ~] = bwlabel(segim1, 4);    
            segims(:, :, k + 1) = uint16(segim1);
            
        end

        for k = 1:size(ims, 3) - 1
            % match_frames contains the mapping algorithm between two
            % frames
                          
            [cell_roster0, back_vals0, mapping, cost_matrices, intensities] = match_frames(cell_roster0, back_vals0, ims(:, :, k), ims(:, :, k + 1), segims(:, :, k), segims(:, :, k + 1), shifts(k, :)); 
            
            % write data from the newly re-segmented image to file
            imagefilename = inputfilenames{1}{i + k};
            outfilename = ['segmented' filesep imagefilename(1:end - 4) '_regionprops.mat'];
            save(outfilename, 'cell_roster0');    
            
            sparse_cost_matrices = sparse(cost_matrices);
            outfilename = ['matrices' filesep imagefilename(1:end - 4) '_costmatrices.mat'];
            save(outfilename, 'sparse_cost_matrices');
            
            outfilename = ['matrices' filesep imagefilename(1:end - 4) '_intensities.mat'];
            save(outfilename, 'intensities');
            
            % update lineage matrix
            lineages = update_lineages_quick(lineages, reduce_match_matrix(mapping), i + k - 1);  
        end         
        
        carryover_im = ims(:, :, end);
        carryover_segim = segims(:, :, end);
        iteration = iteration + 1;
        disp([dir_name ':Segmented ' num2str(processnb*(iteration - 1) + 1) '/' num2str(length(imagefilenames)) ' images']);
    end    
    
    % segment the last par_op images
    par_op = mod(length(imagefilenames) - 1, processnb);
    i = processnb*(iteration - 1) + 1;
    
    % read in 'par_op' images and segment them in parallel
    ims = zeros(imx, imy, par_op + 1, 'uint16');
    segims = ims;
    
    ims(:, :, 1) = carryover_im; segims(:, :, 1) = carryover_segim;
    
    for k = 1:par_op 
        for j = 1:numchannels
            imagefilenames = inputfilenames{1, j};
            ims(:, :, k + 1) = ims(:, :, k + 1) + imread(imagefilenames{i + k});
        end
        [segim1, ~] = bwlabel(laplaceimagesegment_weak(double(ims(:, :, k + 1))./min_filt), 4);
        segims(:, :, k + 1) = uint16(segim1);
    end
    
    % calculates relative shift between two images
    shifts = zeros(par_op, 2);

    for k = 1:size(ims, 3) - 1
        % match_frames contains the mapping algorithm between two
        % frames
        [cell_roster0, back_vals0, mapping, cost_matrices, intensities] = match_frames(cell_roster0, back_vals0, ims(:, :, k), ims(:, :, k + 1), segims(:, :, k), segims(:, :, k + 1), shifts(k, :));
        
        % write data from the newly re-segmented image to file
        imagefilename = inputfilenames{1}{i + k};
        outfilename = ['segmented' filesep imagefilename(1:end - 4) '_regionprops.mat'];
        save(outfilename, 'cell_roster0');
        
        sparse_cost_matrices = sparse(cost_matrices);
        outfilename = ['matrices' filesep imagefilename(1:end - 4) '_costmatrices.mat'];
        save(outfilename, 'sparse_cost_matrices');
        
        outfilename = ['matrices' filesep imagefilename(1:end - 4) '_intensities.mat'];
        save(outfilename, 'intensities');
        
        % update lineage matrix
        lineages = update_lineages(lineages, reduce_match_matrix(mapping), i + k - 1);
    end
    
    outfilename = ['cellproperties' filesep varargin{1} '_lineages.mat'];
    save(outfilename, 'lineages');
    disp([dir_name ':Segmented ' num2str(length(imagefilenames)) '/' num2str(length(imagefilenames)) ' images']);
end


function shifts = find_xcorr(segim0, ~)
    [ind1, ind2, ~] = find(segim0);
    epicenter = [ind1(length(ind1)/2), ind2(length(ind2)/2)];
end

 
%% SEGMENTATION CORE FUNCTION
function seg_im = laplaceimagesegment_weak(back_norm_im)

    global imx thres_offset
    bin = 1024/imx; blur_size = 21; se = strel('square', 3);
    
    % segmentation, roughly adapted from steve eddins' online version
    g_filt = fspecial('gaussian', floor(blur_size/bin), 0.3*(blur_size/(2*bin) - 1) + 0.8);
    gauss_im = imfilter(back_norm_im, g_filt, 'replicate', 'conv');
    max_im_gauss = imregionalmax(gauss_im);
            
    perimeter = zeros(1, 150);
    for k = 1:150
        perimeter(k) = length(find(edge(gauss_im, 'canny', 0.3 - 0.002*(k - 1))));
    end
    [~, ind] = max(smooth(diff(perimeter)));
    
    
%     f1 = figure();
%     plot(0.3 - 0.005*((1:60) - 1), perimeter, 'o')
%     close(f1);
    
%     canny = edge(gauss_im, 'canny', 0.3 - 0.002*(ind - 2));
    canny = edge(gauss_im, 'canny', 0.01);
    cannyplus = double(canny);
    cannyplus(canny == 0) = max(max(back_norm_im)) + 1;
    cannyplus(canny == 1) = back_norm_im(canny == 1);
    canny_filt = ordfilt2(double(cannyplus), 1, ones(10), 'symmetric');
    
    canny_sub = canny_filt;
    canny_sub(canny_filt == max(max(cannyplus))) = 0;
    canny_sub(canny_sub > 0) = 1;
    canny_sub = imfill(canny_sub, 'holes');
        
%     inv_im = imcomplement(gauss_im);
    inv_im = imcomplement(back_norm_im);
    if max(max(max_im_gauss.*canny_sub)) > 0
        inv_im_const = imimposemin(inv_im, max_im_gauss.*canny_sub);        
        watershed_im = watershed(inv_im_const);
        
        %change threshold
        thres = median(inv_im_const(canny_sub == 0)) - thres_offset;
        inv_im_const(inv_im_const > thres) = 0;
        inv_im_const(inv_im_const < thres) = 1;
        seg_im = inv_im_const.*logical(watershed_im);
        seg_im = seg_im.*canny_sub;
        seg_im = imopen(seg_im, se);
%       seg_im = imopen(imclose(imfill(seg_im, 'holes'), se), se);
        
        seg_im = bwlabel(seg_im, 4);
        props = regionprops(seg_im);
        areas = [props.Area];
        for k = 1:length(areas)
            if areas(k) > 1000
                seg_im(seg_im == k) = 0;
            end
        end
    else
        seg_im = 0*inv_im;
        seg_im(1:2, 1:2) = 1;
    end
    seg_im = logical(seg_im);
    
 end


function [cell_roster1, back_vals1, combined_mapping, cost_matrices, intensities] = match_frames(cell_roster0, back_vals0, im0, im1, segim0, segim1, shift)

    global dt it imx dt_max
    bin = 1024/imx;
    nbcell0 = length(cell_roster0);
    
    if nbcell0 > 200
        dt = dt_max/sqrt(nbcell0/100);
    else
        dt = dt_max;
    end
    
    % label second image
    nbcell1 = max(max(segim1));
    cell_roster1 = regionprops(segim1, im1, 'Area', 'BoundingBox', 'Centroid', 'Image', 'Orientation', 'Eccentricity', 'WeightedCentroid', 'MajorAxisLength');
    back_vals1 = zeros(length(cell_roster1), 1); centroids = {cell_roster1.Centroid};
    for k = 1:length(cell_roster1)
        back_vals1(k) = find_background(centroids{k}, im1, segim1, bin);
    end
    
    % calculate background-subtracted-intensities for all cells
    intensity0 = zeros(nbcell0, 1);
    for k = 1:nbcell0
        intensity0(k) = calc_intensity(im0, cell_roster0, back_vals0, k);
    end
        
    intensity1 = zeros(nbcell1, 1);
    for k = 1:nbcell1
        intensity1(k) = calc_intensity(im1, cell_roster1, back_vals1, k);
    end    
           
    intensities = {intensity0, intensity1};
    
    % calculate overlap, distances, difference in intensity between objects in the two labeled images
    [overlap_matrix, distance_matrix, intensity_matrix] = calc_costmatrices(cell_roster0, cell_roster1, segim0, segim1, intensity0, intensity1, shift);
    combined_mapping = munkres_mapping_with_divisions(overlap_matrix, distance_matrix, intensity_matrix, intensity0, intensity1, dt, it);
    cost_matrices = [combined_mapping;overlap_matrix; distance_matrix; intensity_matrix];
end

function [overlapmatrix, distancematrix, intensitymatrix] = calc_costmatrices(cell_roster0, cell_roster1, segim0, segim1, intensity0, intensity1, shift)

    global dt

    nbcell0 = length(cell_roster0); nbcell1 = length(cell_roster1);
    overlapmatrix = zeros(nbcell0, nbcell1);
    distancematrix = overlapmatrix + dt;
    intensitymatrix = overlapmatrix;
    delrow = shift(1); delcol = shift(2);
    for k = 1:nbcell0
        centroid0 = cell_roster0(k).Centroid;
        centroids1 = [cell_roster1.Centroid]; 
        centroids1 = transpose(reshape(centroids1, 2, length(centroids1)/2));
        centroids1 = centroids1 - [delcol*ones(size(centroids1, 1), 1) delrow*ones(size(centroids1, 1), 1)];
        distancematrix(k, :) = transpose(sqrt((centroid0(1) - centroids1(:, 1)).^2 + (centroid0(2) - centroids1(:, 2)).^2));
        intensitymatrix(k, :) = 2*abs(intensity0(k) - intensity1)./(intensity0(k) + intensity1);

        [rows, cols] = find(segim0 == k);
        if (rows(1) + delrow)*(cols(1) + delcol) > 0
            overlapping_region = segim1(rows + delrow, cols + delcol);
        else
            overlapping_region = 0;
        end
        
        overlapping_cells = unique(nonzeros(overlapping_region));
        for p = 1:length(overlapping_cells)
            overlapmatrix(k, overlapping_cells(p)) = length(nonzeros(overlapping_region == overlapping_cells(p)))/(cell_roster0(k).Area + cell_roster1(overlapping_cells(p)).Area);
        end
    end
    
    distancematrix(distancematrix > dt) = dt;
    
end

function bound = cart2mat(bound)
    temp = bound;
    bound(1) = temp(2);
    bound(2) = temp(1);
end 

function [range_row, range_col] = expand_box(bound0, bound1)
    range_row = max(bound0(1) + bound0(3) - 1, bound1(1) + bound1(3) - 1) - min(bound0(1), bound1(1)) + 1;
    range_col = max(bound0(2) + bound0(4) - 1, bound1(2) + bound1(4) - 1) - min(bound0(2), bound1(2)) + 1;
end


function back_sub_intensity = calc_intensity(im0, cell_roster0, back_vals0, index)
    cellim0 = double(find_subimage(im0, cell_roster0(index).BoundingBox, cell_roster0(index).Image));
    back_sub_intensity = sum(nonzeros(cellim0)) - length(nonzeros(cellim0))*back_vals0(index);
end

