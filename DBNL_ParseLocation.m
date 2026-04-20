function [hemisphere, lobe, structure] = DBNL_ParseLocation(electrode_loc)
    % Initialize outputs
    hemisphere = NaN;
    lobe = NaN;
    structure = NaN;

    % Split the string into parts using '-' as a delimiter
    parts = split(electrode_loc, '-');
    
    % Check if the string has exactly 2 parts
    if numel(parts) == 2
        hemisphereCode = parts{1};
        structure = parts{2};

        % Parse the hemisphere
        switch hemisphereCode
            case {'LT', 'LP', 'LO', 'LF', 'L'}
                hemisphere = 'Left';
            case {'RT', 'RP', 'RO', 'RF', 'R'}
                hemisphere = 'Right';
            case 'Z'
                hemisphere = 'None'; % Unknown hemisphere
        end

        % Parse the lobe
        if ~ismember(hemisphereCode, {'Z'}) % Exclude non-brain locations
            switch hemisphereCode
                case {'LT', 'RT'}
                    lobe = 'Temporal';
                case {'LF', 'RF'}
                    lobe = 'Frontal';
                case {'LP', 'RP'}
                    lobe = 'Parietal';
                case {'LO', 'RO'}
                    lobe = 'Occipital';
                otherwise
                    lobe = 'None'; % For structures without a lobe
            end
        end

        % Special case: Z-Unknown
        if strcmp(hemisphereCode, 'Z') || strcmp(structure, 'Unknown')
            hemisphere = NaN;
            lobe = NaN;
            structure = NaN;
        end
    end
end