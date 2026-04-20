function post_data = DBNL_preprocess(data, fs, ref, filt, removestim, TTL, suppressidx, LR_idx)
    % DBNL_preprocess Preprocesses the input data by removing stimulation artifacts, filtering, and re-referencing.
    %
    % INPUTS:
    % data         - Matrix of raw data (channels x timepoints)
    % fs           - Sampling frequency
    % ref          - Reference scheme ('CAR', 'BP', 'LR', or none)
    % filt         - Filter parameters (empty if no filtering)
    % removestim   - Flag indicating whether to remove stimulation artifacts (1 = yes, 0 = no)
    % TTL          - TTL signal for identifying stimulation artifacts
    % suppressidx  - Indices of channels to suppress during re-referencing
    % LR_idx       - Optional input that specifies an index to locally reference the data to
    %
    % OUTPUT:
    % post_data    - Preprocessed data matrix

    if nargin < 8
        LR_idx = [];
    end

    % Remove stimulation artifacts if specified
    if removestim == 1
        TTL_timepoints = find(TTL == 1);
        if size(TTL_timepoints, 2) == 1
            % Transpose the row vector to a column vector (N x 1)
            TTL_timepoints = TTL_timepoints';
        end
        uni_TTL = TTL_timepoints([diff(TTL_timepoints) ~= 1, true]);
        TTL_before = round(15 / 1000 * fs);
        TTL_after = round(15 / 1000 * fs);
        for TTLnum = 1:numel(uni_TTL)
            Before_data = data(:, uni_TTL(TTLnum) - TTL_before : uni_TTL(TTLnum));
            taper = linspace(0, 1, size(Before_data, 2));
            Before_taper = Before_data .* flip(taper);
            After_data = data(:, uni_TTL(TTLnum) + TTL_after : uni_TTL(TTLnum) + TTL_after*2);
            After_taper = After_data .* taper;
            data_replace = Before_taper + After_taper;
            data(:, uni_TTL(TTLnum) : uni_TTL(TTLnum) + TTL_after) = data_replace;
        end
    end

    % Remove 60Hz line noise
    noline_data = DBNL_remove60hz(data, fs);
    
    % Preallocate filtered data matrix
    [num_channels, num_timepoints] = size(data);
    filt_data = nan(num_channels, num_timepoints);
    
    % Filter the data if filter parameters are provided
    parfor chan_num = 1:num_channels
        if isempty(filt)
            filt_data(chan_num, :) = noline_data(chan_num, :);
        else
            [filt_data(chan_num, :), ~, ~] = Butterworth_Hilbert_LR(noline_data(chan_num, :), fs, filt);
        end
    end

    % Re-reference the data based on the specified scheme
    switch ref
        case 'CAR'
            post_data = nan(size(data));
            temp_data = filt_data;
            temp_data(suppressidx, :) = NaN;
            mean_temp_data = mean(temp_data, 1, "omitnan");
            for timepoint = 1:num_timepoints
                post_data(:, timepoint) = filt_data(:, timepoint) - mean_temp_data(timepoint);
            end
            post_data(suppressidx,:) = NaN;
        case 'CARred'
            post_data = nan(size(data));
            temp_data = filt_data(LR_idx,:);
            mean_temp_data = mean(temp_data, 1, "omitnan");
            for timepoint = 1:num_timepoints
                post_data(:, timepoint) = filt_data(:, timepoint) - mean_temp_data(timepoint);
            end
            post_data(suppressidx,:) = NaN;
        case 'BP'
            post_data = nan(size(data));
            for bipol_num = 1:num_channels-1
                post_data(bipol_num, :) = filt_data(bipol_num, :) - filt_data(bipol_num + 1, :);
            end
            post_data(suppressidx, :) = NaN;
        case 'BBP'
            post_data = nan(num_channels,num_channels,num_timepoints);
            filt_data(suppressidx,:) = NaN;
            for bipol_num = 1:num_channels
                for bipol_num2 = 1:num_channels
                    post_data(bipol_num, bipol_num2, :) = filt_data(bipol_num, :) - filt_data(bipol_num2, :);
                end
            end
        case 'LP'
            post_data = nan(size(data));
            for laplace_num = 2:num_channels-1
                post_data(laplace_num, :) = filt_data(laplace_num, :) - (filt_data(laplace_num - 1, :) + filt_data(laplace_num + 1, :))./2;
            end
            post_data(suppressidx, :) = NaN;
        case 'LR'
            post_data = filt_data - filt_data(LR_idx, :);
            post_data(suppressidx, :) = NaN;
        otherwise
            filt_data(suppressidx, :) = NaN;
            post_data = filt_data;
    end
    post_data(num_channels+1:end,:) = [];
end