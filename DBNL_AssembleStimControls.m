function [Stim_Controls, Stim_Control_Labels] = DBNL_AssembleStimControls(stim_conditions, stimidx)
    % Split the stim_conditions string into individual components
    stim_order = split(stim_conditions, ', ');

    % Initialize variables for storing the extracted information
    periodicity = strings(size(stim_order));
    current = strings(size(stim_order));

    % Loop through each component and parse both periodicity and current
    for i = 1:length(stim_order)
        if contains(stim_order(i), 'NP')
            periodicity(i) = "Non-Periodic";
        elseif contains(stim_order(i), 'P')
            periodicity(i) = "Periodic";
        end

        if contains(stim_order(i), 'HA')
            current(i) = "High";
        elseif contains(stim_order(i), 'LA')
            current(i) = "Low";
        end
    end

    % Check for high amplitude condition before the given index
    if stimidx > 1 % Ensure stimidx is valid
        Prior = any(current(1:stimidx-1) == "High");
    else
        Prior = false; % No previous conditions to check
    end

    % Create output controls and labels
    Stim_Controls = {periodicity(stimidx), current(stimidx), Prior, stimidx};
    Stim_Control_Labels = {'Pattern', 'Current', 'HA_Prior', 'Stim_Order'};
end