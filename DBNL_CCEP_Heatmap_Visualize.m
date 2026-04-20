function [sig_heatmap, amp_heatmap] = DBNL_CCEP_Heatmap_Visualize(CCEP_mat1, CCEP_mat2, timepoints, CRP_analysis_idx, analysis_chan, FS_Label, suppress_idx, SR, plot_title)
    % DBNL_CCEP_Heatmap_Visualize Visualizes the significant differences in CCEP data as a heatmap.
    %
    % INPUTS:
    % CCEP_mat1         - 3D matrix of CCEP data for condition 1 (channels x timepoints x trials)
    % CCEP_mat2         - 3D matrix of CCEP data for condition 2 (channels x timepoints x trials)
    % timepoints        - Vector of time points corresponding to the CCEP data
    % CRP_analysis_idx  - Index of the starting timepoint for CRP analysis
    % analysis_chan     - Vector of indices indicating which channels to analyze
    % FS_Label          - Cell array of labels for the functional sites
    % suppress_idx      - Indices of channels to suppress in the significance analysis
    % plot_title        - String defining the title for the figure e.g. "NP HA", "NP LA", "P HA", "P LA"
    %
    % OUTPUT:
    % sig_heatmap       - Heatmap matrix showing significant differences between conditions
    %
    % USE CASE:
    % This function is used to visualize significant differences in CCEP data between two condi
    % by generating a heatmap. The heatmap highlights significant time points and channels bad on
    % a specified threshold.

    % Sort labels and indices for functional sites
    %[sort_label, sort_idx] = sortrows(FS_Label');
    
    % Find the subsampling index stopping at 500ms
    subsample_timepoints = 0.015:0.001:0.2;
    for i = 1:numel(subsample_timepoints)
        [~, timepoints_sub_idx(i)] = min(abs(timepoints - subsample_timepoints(i)));
    end
    
    % Select timepoints for analysis, sampled at every ms interval
    %timepoints_sub_idx = CRP_analysis_idx:round(SR/1000):stop_idx;
    p_thresh = 0.05 / numel(timepoints_sub_idx);
    %p_thresh = 0.01;

    % Compare timepoints between the two conditions
    [p_timepoints] = DBNL_CCEP_timepoint_blockcompare(CCEP_mat1, CCEP_mat2, timepoints_sub_idx);
    
    % Suppress specified indices in the significance analysis
    p_timepoints(suppress_idx, :) = 1;

    % Determine significant time points based on the threshold
    p_timepoints_sig = zeros(size(p_timepoints));
    p_timepoints_sig(p_timepoints < p_thresh) = 1;

    % Calculate the average CCEP data for both conditions
    CCEP_mat1_avg = mean(CCEP_mat1, 3);
    CCEP_mat2_avg = mean(CCEP_mat2, 3);

    % Compute the difference between the two conditions
    CCEP_Diff = (CCEP_mat2_avg(:, timepoints_sub_idx) - CCEP_mat1_avg(:, timepoints_sub_idx));
    CCEP_Diff(suppress_idx,:,:) = 0;

    [CCEP_pre_amp, ~] = envelope_skip_allnan(CCEP_mat1_avg(:,timepoints_sub_idx)',300,'rms');
    [CCEP_post_amp, ~] = envelope_skip_allnan(CCEP_mat2_avg(:,timepoints_sub_idx)',300,'rms');
    CCEP_amp_Diff = CCEP_post_amp - CCEP_pre_amp;
    CCEP_amp_Diff = CCEP_amp_Diff';
    CCEP_amp_Diff(suppress_idx,:,:) = 0;
    amp_heatmap = CCEP_amp_Diff.*p_timepoints_sig;

    % Correct heatmap labels and indices
    [heatmap_label, heatmap_idx] = DBNL_correctheatmap(sort_label);
    [~, ia] = unique(sort_label);

    % Generate the significance heatmap
    figure();
    sig_heatmap = p_timepoints_sig .* CCEP_Diff;
    sig_heatmap(isnan(sig_heatmap)) = 0;
    %sig_heatmap = CCEP_Diff;
    imagesc(timepoints(timepoints_sub_idx) * 1000, 1:numel(analysis_chan), sig_heatmap(sort_idx, :));
    yticks(heatmap_idx);
    yticklabels(heatmap_label);
    %clim([-150 150]);
    colormap(redblue);
    colorbar();
    title(plot_title)
    
    % Add vertical lines to separate different groups
    hline(ia(2:end) - 0.5, 'k-');
    hold on;
    
    % Add labels and colorbar
    xlabel('Time From Stimulation (ms)');
    cb = colorbar();
    ylabel(cb, 'Post-CCEP - Pre-CCEP (uV)', 'FontSize', 16, 'Rotation', 270);
end