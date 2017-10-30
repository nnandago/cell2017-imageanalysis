function processfluor_MAC(dir_name, segchannel, bin, varargin)
% processes images for fluorescence extraction for masks in segchannel
% can be called with extra inputs in the form {fluorescence_channel,
% background_channel}, where maskes in 'background_channel' are used for
% background correction.

    global processnb imx imy

    % number of processes to run in parallel
    processnb = 2; %str2double(getenv('NUMBER_OF_PROCESSORS'));
    imx = 1024/bin; imy = 1024/bin;
        
    optargin = size(varargin, 2);
    num_channels = 0;
    if optargin >= 1
        num_channels = optargin;
    end
    
     
    if exist(dir_name, 'dir')
        cd(dir_name);
    else
        error('Invalid directory. Please try again');
    end
    
    if ~exist('segmented', 'dir')
        error('Please run segtrack.m first');
    end
    
    vars = load(['cellproperties' filesep segchannel, '_lineages.mat']);
    lineages = vars.lineages + 1;
    
    for j = 1:num_channels + 1
        if j == 1
            channel = segchannel;
            background_channel = segchannel;
            channel_diff = 0;
        else
            arg = varargin{j - 1};
            channel = char(arg{1});
            background_channel = char(arg{2});
            if strcmp(channel, background_channel) == 1
                channel_diff = 0;
            else
                channel_diff = 1;
            end
        end
        
        imagefilenames = dir(['*' channel '*.TIF']);
        propfilenames = dir(['segmented' filesep '*', segchannel, '*regionprops.mat']);   
        
        indices2 = zeros(length(propfilenames), 1);
        for k = 1:length(propfilenames)
            fn = propfilenames(k).name;
            fn = fn(regexp(fn, segchannel):end);
            start_ind = regexp(fn, '_t');
            indices2(k) = str2double(fn(start_ind + 2:end - 16));
        end
        
        [~, sort_order] = sort(indices2);
        propfilenames = {propfilenames(sort_order).name};
        if j == 1
            maxfiles = length(propfilenames);
        end
        
        % use masks in background_channel if channel is different from
        % background_channel        
        if channel_diff
            propfilenames_back = dir(['segmented' filesep '*', background_channel, '*regionprops.mat']);
            propfilenames_back = {propfilenames_back(sort_order).name};
        else
            propfilenames_back = propfilenames;
        end
        
        indices1 = zeros(length(imagefilenames), 1);
        for k = 1:length(imagefilenames)
            fn = imagefilenames(k).name;
            fn = fn(regexp(fn, channel):end);
            start_ind = regexp(fn, '_t');
            indices1(k) = str2double(fn(start_ind + 2:end - 4));
        end
            
        [~, sort_order] = sort(indices1);
        imagefilenames = {imagefilenames(sort_order).name};
        imagefilenames = imagefilenames(1:maxfiles);
        
                                     
        cell_indices = cell(size(lineages, 2), 1);
        for k = 1:size(lineages, 2)
            cell_indices{k} = find(lineages(:, k));
        end
                
        
        iteration = 1; tot_iterations = floor((maxfiles - 1)/processnb);
%         h = waitbar(0, ['Channel: ' channel '. Processed: ' num2str(iteration) '/' num2str(maxfiles) ' image(s)']);
        
        tot_fluorescence = zeros(size(lineages));
        mean_background = zeros(size(lineages));
        areas = zeros(size(lineages));
        med_fluorescence = zeros(size(lineages));
        
        while iteration <= tot_iterations
            
            i = processnb*(iteration - 1);
                       
            for k = 1:processnb
                roster = load(['segmented' filesep propfilenames{i + k}]);
                back_roster = load(['segmented' filesep propfilenames_back{i + k}]);
                im = imread(imagefilenames{i + k});
                properties = extractfluor(roster, back_roster, im, bin);
                tot_fluorescence(:, i + k) = properties(lineages(:, i + k), 1);
                mean_background(:, i + k) = properties(lineages(:, i + k), 2);
                med_fluorescence(:, i + k) = properties(lineages(:, i + k), 3);
                areas(:, i + k) = properties(lineages(:, i + k), 4);                               
            end   
            
            iteration = iteration + 1;
%             waitbar((processnb*(iteration - 1) + 1)/length(imagefilenames), h, ['Channel: ' channel '. Processed: ' num2str(processnb*(iteration - 1) + 1) '/' num2str(length(imagefilenames)) ' images'])
        end 
        
        par_op = mod(length(imagefilenames) - 1, processnb) + 1;
        i = processnb*(iteration - 1);
        
        for k = 1:par_op %, par_op)
            roster = load(['segmented' filesep propfilenames{i + k}]);
            back_roster = load(['segmented' filesep propfilenames_back{i + k}]);
            im = imread(imagefilenames{i + k});
            properties = extractfluor(roster, back_roster, im, bin);
            tot_fluorescence(:, i + k) = properties(lineages(:, i + k), 1);
            mean_background(:, i + k) = properties(lineages(:, i + k), 2);
            med_fluorescence(:, i + k) = properties(lineages(:, i + k), 3);
            areas(:, i + k) = properties(lineages(:, i + k), 4);
        end
        
        save(['cellproperties' filesep channel, '_' segchannel '_data.mat'], 'tot_fluorescence', 'mean_background', 'med_fluorescence', 'areas')  
%         waitbar(1, h, ['Channel: ' channel '. Processed ' num2str(length(imagefilenames)) '/' num2str(length(imagefilenames)) ' images'])
%         close(h)
    end
end

function props = extractfluor(roster_struct, mask_roster_struct, im, bin)
    
    global imx imy

    mask_cell_roster = mask_roster_struct.cell_roster0;
    cell_roster = roster_struct.cell_roster0;
    props = zeros(length(cell_roster) + 1, 4);
    segim = zeros(imx, imy);
    
    for k = 1:length(mask_cell_roster)
        bounds = uint16(mask_cell_roster(k).BoundingBox);
        cellim = mask_cell_roster(k).Image;
        [upperbound_x, upperbound_y, diff_x, diff_y] = find_bounds([bounds(2), bounds(1)], size(cellim), [imx imy]);
        segim(bounds(2):upperbound_x, bounds(1):upperbound_y) = segim(bounds(2):upperbound_x, bounds(1):upperbound_y) + cellim(1:size(cellim, 1) - diff_x, 1:size(cellim, 2) - diff_y);
    end
        
    for k = 2:size(props, 1)
        fluor_subim = find_subimage(im, cell_roster(k-1).BoundingBox, cell_roster(k-1).Image);
        props(k, 1) = sum(sum(fluor_subim));
        props(k, 2) = find_background(cell_roster(k-1).Centroid, im, segim, bin);  
        props(k, 3) = median(median(double(fluor_subim), 1));
        props(k, 4) = cell_roster(k-1).Area;
    end
end