function [CCEP_mat, timepoints] = DBNL_CCEP_segmentation(filt_data, TTL, bounds, fs, suppress_idx, ref, exclude_trials)
    % CCEP_segmentation Segments filtered data based on TTL signals and specified bounds.
    %This file assumes the download of the CARLA function "ccep_CARVariance.m"
    %
    % INPUTS:
    % filt_data     - Filtered data matrix (channels x timepoints)
    % TTL           - TTL signal for identifying stimulation artifacts
    % bounds        - Time bounds for segmentation (e.g., [-0.5, 2])
    % fs            - Sampling frequency
    % suppress_idx  - Indices of channels to suppress during processing
    % useCARLA      - Flag to indicate whether to use CARLA for variance reduction (1 = yes, 0 = no)
    % exclude_trials- Vector of trial indices to exclude from analysis
    %
    % OUTPUTS:
    % CCEP_mat      - Segmented CCEP data (channels x timepoints x responses)
    % timepoints    - Time points vector for the segmented data

    % Find unique TTL time points (start of each stimulation)
    TTL_timepoints = find(TTL == 1);
    if size(TTL_timepoints, 2) == 1
            % Transpose the row vector to a column vector (N x 1)
            TTL_timepoints = TTL_timepoints';
    end
    uni_TTL = TTL_timepoints([diff(TTL_timepoints) ~= 1, true]);
    
    % Calculate bounds in samples
    low_bound = round(fs * bounds(1));
    up_bound = round(fs * bounds(2));
    base_bound_up = round(fs*0.015);
    base_bound_low = round(fs*0.2);
    
    % Generate timepoints vector
    timepoints = bounds(1):1/fs:bounds(2);
    
    % Preallocate CCEP_mat for efficiency
    num_responses = numel(uni_TTL) - 1 - numel(exclude_trials);
    num_channels = size(filt_data, 1);
    num_timepoints = numel(timepoints);

    if ndims(filt_data) == 2
        CCEP_mat = zeros(num_channels, num_timepoints, num_responses);
    
        % Initialize response index
        resp_idx = 1;
        for CCEP_resp = 1:numel(uni_TTL)
            if ismember(CCEP_resp, exclude_trials)
                continue; % Skip excluded trials
            end
    
            norm = 0;
            if bounds(1) < 0
                % Normalize based on pre-stimulation period
                %norm = 0;
                norm = mean(filt_data(:, uni_TTL(CCEP_resp) - base_bound_low : uni_TTL(CCEP_resp)-base_bound_up), 2);
            end
            % Segment and normalize the data
            CCEP_mat(:, :, resp_idx) = filt_data(:, uni_TTL(CCEP_resp) + low_bound : uni_TTL(CCEP_resp) + up_bound) - norm;
            resp_idx = resp_idx + 1;
        end
    
        % Apply CARLA for variance reduction if specified
        if ref == "CARLA"
            [CCEP_mat, ~] = ccep_CARVariance(timepoints, CCEP_mat, fs, suppress_idx);
        end
    elseif ndims(filt_data) == 3
        CCEP_mat = zeros(num_channels, num_channels, num_timepoints, num_responses);
    
        % Initialize response index
        resp_idx = 1;
        for CCEP_resp = 1:numel(uni_TTL)
            if ismember(CCEP_resp, exclude_trials)
                continue; % Skip excluded trials
            end
    
            norm = 0;
            if bounds(1) < 0
                % Normalize based on pre-stimulation period
                norm = mean(filt_data(:, :, uni_TTL(CCEP_resp) - base_bound_low : uni_TTL(CCEP_resp)-base_bound_up), 3);
            end
            % Segment and normalize the data
            CCEP_mat(:, :, :, resp_idx) = filt_data(:, :, uni_TTL(CCEP_resp) + low_bound : uni_TTL(CCEP_resp) + up_bound) - norm;
            resp_idx = resp_idx + 1;
        end
    
        % Apply CARLA for variance reduction if specified
        if ref == "CARLA"
            [CCEP_mat, ~] = ccep_CARVariance(timepoints, CCEP_mat, fs, suppress_idx);
        end
    end
end