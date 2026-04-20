function [suppress_idx, suppress_idx_no_stim, suppress_idx_no_sub] = DBNL_define_supress_chan(channel_names, search_strings_file)
    % DBNL_define_supress_chan Searches for strings in channel_names and returns the index numbers.
    %
    % INPUTS:
    % channel_names       - String array of channel names (e.g., ["LANT1", "RAOF5", "LSG7", ...])
    % search_strings_file - Filename of the text file containing channel names to exclude
    %
    % OUTPUTS:
    % suppress_idx            - Vector of indices in channel_names that match any of the search strings
    % suppress_idx_no_stim    - Vector of indices in channel_names that match any of the search strings excluding those with keyword "stim"
    % suppress_idx_no_sub    - Vector of indices in channel_names that match any of the search strings excluding those with keyword "sub"

    % Validate inputs
    if ~isstring(channel_names) || ~ischar(search_strings_file)
        error('channel_names must be a string array and search_strings_file must be a character array (string).');
    end

    % Read the search strings from the text file
    data = readtable(search_strings_file, 'Delimiter', ',', 'ReadVariableNames', false, 'Format', '%s%s');

    % Extract search strings and descriptions
    search_strings = string(data.Var1);
    descriptions = string(data.Var2);

    suppress_idx = [];

    for i = 1:length(search_strings)
        pattern = strcat('(^|\W)', search_strings{i}, '($|\W)');
        matches = ~cellfun('isempty', regexp(channel_names, pattern, 'once'));
        suppress_idx = [suppress_idx; find(matches)];
    end
    
    suppress_idx = unique(suppress_idx);

    % Find indices excluding those with keyword "stim"
    no_stim_mask = ~contains(descriptions, 'stim', 'IgnoreCase', true);
    search_strings_no_stim = search_strings(no_stim_mask);
    suppress_idx_no_stim = find(ismember(channel_names, search_strings_no_stim));

    % Find indices excluding those with keyword "sub"
    no_sub_mask = ~contains(descriptions, 'sub', 'IgnoreCase', true);
    search_strings_no_sub = search_strings(no_sub_mask);
    suppress_idx_no_sub = find(ismember(channel_names, search_strings_no_sub));
end