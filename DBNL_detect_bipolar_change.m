function jumpidx = DBNL_detect_bipolar_change(Channel_Names)
    % Initialize an empty array to hold the indices of changes
    jumpidx = [];

    % Loop through the string array
    for i = 1:length(Channel_Names)-2
        % Extract the suffix numeric part of the current and next elements
        currentSuffix = regexp(Channel_Names{i}, '\d+$', 'match');
        nextSuffix = regexp(Channel_Names{i+1}, '\d+$', 'match');
        
        % Convert the suffix to numbers
        currentNumber = str2double(currentSuffix{1});
        nextNumber = str2double(nextSuffix{1});
        
        % Check if the sequence number changes
        if nextNumber ~= currentNumber + 1
            jumpidx = [jumpidx; i+1];
        end
    end
    jumpidx = [jumpidx; length(Channel_Names)];
end