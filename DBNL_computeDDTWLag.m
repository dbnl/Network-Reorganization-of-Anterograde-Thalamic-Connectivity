function [t_eval, lag_periodic_all, lag_nonperiodic_all] = DBNL_computeDDTWLag(CCEP_mat,Order, t)
% computeDDTWLag - Compute instantaneous lag using Derivative DTW (DDTW)
% while skipping channels with NaN data in the evaluation interval.
%
% Syntax:
%   [t_eval, lag_periodic, lag_nonperiodic] = computeDDTWLag(data, t)
%
% Inputs:
%   data - A C x T x 4 matrix where:
%          C = number of channels,
%          T = number of time points,
%          The 4 stimulation conditions are in the order:
%             1: Pre-periodic, 2: Post-periodic,
%             3: Pre-nonperiodic, 4: Post-nonperiodic.
%   t    - A 1 x T time vector.
%
% Outputs:
%   t_eval          - Time vector restricted to the evaluation interval [15, 500] ms.
%   lag_periodic    - Average instantaneous lag (in ms) for the periodic condition.
%   lag_nonperiodic - Average instantaneous lag (in ms) for the nonperiodic condition.
%

eval_idx = (t >= 15) & (t <= 500);
t_eval = t(eval_idx);
dt = t_eval(2) - t_eval(1);
L = length(t_eval);
data = [];
for i = 1:numel(Order)
    data(:,:,i) = mean(CCEP_mat{Order(i)},3);
end

numChannels = size(data, 1);

% Preallocate matrices to store instantaneous lag for each channel.
lag_periodic_all = nan(numChannels, L);
lag_nonperiodic_all = nan(numChannels, L);

for ch = 1:numChannels
    pre_per = squeeze(data(ch, eval_idx, 1));   % Pre-periodic
    post_per = squeeze(data(ch, eval_idx, 2));    % Post-periodic
    
    % If either signal has NaNs in the evaluation interval, skip this channel.
    if any(isnan(pre_per)) || any(isnan(post_per))
        lag_periodic_all(ch, :) = nan(1, L);
    else
        % Compute derivatives.
        dPre_per = gradient(pre_per, dt);
        dPost_per = gradient(post_per, dt);
        
        % Compute DTW on the derivative signals.
        [~, ix, iy] = dtw(dPre_per, dPost_per);
        
        % Remove duplicate indices in ix by averaging corresponding iy values.
        [unique_ix, ~, ic] = unique(ix);
        unique_iy = accumarray(ic, iy, [], @mean);
        
        % Create an interpolant for the warping function.
        warpingFunc = griddedInterpolant(unique_ix, unique_iy, 'linear', 'nearest');
        
        % For each index in the evaluation window, estimate the corresponding index.
        f_est = warpingFunc(1:L);
        
        % Compute instantaneous lag (in indices).
        lag_periodic_all(ch, :) = ((1:L) - f_est) * dt;
    end
    
    pre_nonper = squeeze(data(ch, eval_idx, 3));   % Pre-nonperiodic
    post_nonper = squeeze(data(ch, eval_idx, 4));    % Post-nonperiodic
    
    if any(isnan(pre_nonper)) || any(isnan(post_nonper))
        lag_nonperiodic_all(ch, :) = nan(1, L);
    else
        dPre_nonper = gradient(pre_nonper, dt);
        dPost_nonper = gradient(post_nonper, dt);
        
        [~, ix, iy] = dtw(dPre_nonper, dPost_nonper);
        [unique_ix, ~, ic] = unique(ix);
        unique_iy = accumarray(ic, iy, [], @mean);
        warpingFunc = griddedInterpolant(unique_ix, unique_iy, 'linear', 'nearest');
        f_est = warpingFunc(1:L);
        lag_nonperiodic_all(ch, :) = ((1:L) - f_est) * dt;
    end
end

%% Compute average lag across channels (ignoring NaN rows)
lag_periodic = nanmean(lag_periodic_all, 1);
lag_nonperiodic = nanmean(lag_nonperiodic_all, 1);

end