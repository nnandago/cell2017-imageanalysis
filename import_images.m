function imagefilenames = import_images(dir_name, channel)

        % sort filenames numerically
        imagefilenames = dir([dir_name filesep '*' channel '*TIF']);
        indices = zeros(length(imagefilenames), 1);
        for k = 1:length(imagefilenames)
            fn = imagefilenames(k).name;
            fn = fn(regexp(fn, channel):end);
            start_ind = regexp(fn, '_t');
            indices(k) = str2double(fn(start_ind + 2:end - 4));
        end
        
        [~, sort_order] = sort(indices);
        
        % structure of sorted filenames
        imagefilenames = {imagefilenames(sort_order).name};
end