function [Contact_Controls, Contact_Control_labels] = DBNL_AssembleContactControls(parsecontact, stim_contacts, SOZ, electrodecsv)
    % Extract electrode data
    electrode_loc = electrodecsv.FSLabel_aparc_DKTatlas_aseg;
    electrode_loc = DBNL_autoseg_cleanup(electrode_loc);
    electrode_labels = electrodecsv.Label;

    % Find the location of the parsed contact
    parse_electrode_loc = find(strcmp(electrode_labels, parsecontact));

    % Parse SOZ and stim contacts
    parse_SOZ_loc = [];
    parse_stim_loc = [];
    SOZ = split(SOZ, ', ');
    stim_contacts = split(stim_contacts, ', ');
    for i = 1:numel(SOZ)
        parse_SOZ_loc(i) = find(strcmp(electrode_labels, SOZ(i)));
    end
    for j = 1:numel(stim_contacts)
        parse_stim_loc(j) = find(strcmp(electrode_labels, stim_contacts(j)));
    end

    % Calculate distances
    [c2c_distance] = DBNL_c2clocalize(electrodecsv);
    SOZ_dist = min(c2c_distance(parse_electrode_loc, parse_SOZ_loc));
    stim_dist = min(c2c_distance(parse_electrode_loc, parse_stim_loc));

    % Parse electrode locations
    [hemisphere_cont, lobe_cont, structure_cont] = DBNL_ParseLocation(electrode_loc(parse_electrode_loc));
    [hemisphere_stim, ~, ~] = DBNL_ParseLocation(electrode_loc(parse_stim_loc(1)));

    % Ensure hemisphere_cont is a string or character vector
    if isnumeric(hemisphere_cont)
        hemisphere_cont = string(num2str(hemisphere_cont));
    elseif iscell(hemisphere_cont) && numel(hemisphere_cont) == 1
        hemisphere_cont = string(hemisphere_cont{1});
    elseif ischar(hemisphere_cont)
        hemisphere_cont = string(hemisphere_cont);
    elseif ~isstring(hemisphere_cont)
        error('Unsupported type for hemisphere_cont.');
    end

    % Ensure lobe_cont is a string or character vector
    if isnumeric(lobe_cont)
        lobe_cont = string(num2str(lobe_cont));
    elseif iscell(lobe_cont) && numel(lobe_cont) == 1
        lobe_cont = string(lobe_cont{1});
    elseif ischar(lobe_cont)
        lobe_cont = string(lobe_cont);
    elseif ~isstring(lobe_cont)
        error('Unsupported type for lobe_cont.');
    end

    % Ensure hemisphere_stim is a string or character vector
    if isnumeric(hemisphere_stim)
        hemisphere_stim = string(num2str(hemisphere_stim));
    elseif iscell(hemisphere_stim) && numel(hemisphere_stim) == 1
        hemisphere_stim = string(hemisphere_stim{1});
    elseif ischar(hemisphere_stim)
        hemisphere_stim = string(hemisphere_stim);
    elseif ~isstring(hemisphere_stim)
        error('Unsupported type for hemisphere_stim.');
    end

    % Convert to categorical
    hemisphere_cont_cat = categorical(hemisphere_cont, {'Left', 'Right'});
    lobe_cont_cat = categorical(lobe_cont, {'Temporal', 'Occipital', 'Frontal', 'Parietal', 'None'});
    hemisphere_stim_cat = categorical(hemisphere_stim, {'Left', 'Right'});

    % Return as normalized representations
    Contact_Controls = {hemisphere_cont_cat, lobe_cont_cat, structure_cont, hemisphere_stim_cat, SOZ_dist, stim_dist};
    Contact_Control_labels = {'Cont_hemi', 'Cont_lobe', 'Cont_struct', 'Stim_hemi', 'SOZ_dist', 'stim_dist'};
end