function [p_timestamps] = DBNL_CCEP_timepoint_blockcompare_phase(CCEP_mat1, CCEP_mat2, timepoints_sub_idx)
    % DBNL_CCEP_timepoint_blockcompare Performs statistical comparison of CCEP data at specified time points.
    %
    % INPUTS:
    % CCEP_mat1         - 3D matrix of CCEP data for condition 1 (channels x timepoints x trials)
    % CCEP_mat2         - 3D matrix of CCEP data for condition 2 (channels x timepoints x trials)
    % timepoints_sub_idx - Vector of indices indicating which time points to compare
    %
    % OUTPUT:
    % p_timestamps      - Matrix of p-values for each channel and time point comparison
    %
    % USE CASE:
    % This function performs a two-sample t-test comparing CCEP data between two conditions
    % at specified time points for each channel. It returns a matrix of p-values indicating
    % the statistical significance of the differences.

    % Preallocate the p_timestamps matrix for efficiency
    num_channels = size(CCEP_mat1, 1);
   
    num_timepoints = numel(timepoints_sub_idx);
    p_timestamps = zeros(num_channels, num_timepoints);

    % Perform the two-sample t-test for each channel and time point
    for i = 1:num_channels
        for j = 1:num_timepoints
            % Extract data for the current channel and time point from both conditions
            data1 = squeeze(CCEP_mat1(i, timepoints_sub_idx(j), :));
            data2 = squeeze(CCEP_mat2(i, timepoints_sub_idx(j), :));
            
            % Perform the t-test and store the p-value
            [p, ~, ~] = circ_kuipertest(data1, data2);
            p_timestamps(i, j) = p;
        end
    end
end