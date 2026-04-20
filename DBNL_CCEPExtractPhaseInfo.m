function [CCEP_mat_phase, CCEP_mat_amp] = DBNL_CCEPExtractPhaseInfo(CCEP_mat)
cell_size = numel(CCEP_mat);
chan_num = numel(CCEP_mat{1}(:,1,1));
for count = 1:cell_size
    trial_data = CCEP_mat{count};
    stim_num = numel(trial_data(1,1,:));
    for chan = 1:chan_num
        for stim = 1:stim_num
            ana_signal = squeeze(trial_data(chan,:,stim));
            trans_signal = hilbert(ana_signal);
            phs_fts = angle(trans_signal);
            abs_fts = abs(trans_signal);
            CCEP_mat_phase{count}(chan,:,stim) = phs_fts;
            CCEP_mat_amp{count}(chan,:,stim) = abs_fts;
        end
    end
end
end