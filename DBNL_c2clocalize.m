function [c2c_distance, found_idx] = DBNL_c2clocalize(electrodecsv, targetlist)
    % DBNL_c2clocalize
    % Calculate contact-to-contact distances for all requested targets.
    %
    % If a target in targetlist is not found in electrodecsv.Label, its
    % corresponding row/column in c2c_distance is left as NaN.
    %
    % Inputs:
    %   electrodecsv : table containing Label, Coord_x, Coord_y, Coord_z
    %   targetlist   : cell array or string array of target contact labels
    %
    % Outputs:
    %   c2c_distance : [N x N] distance matrix in same order as targetlist
    %   found_idx    : [N x 1] row index into electrodecsv for each target;
    %                  NaN if target was not found

    contacttable = electrodecsv;
    targetlist = string(targetlist);
    labels = string(contacttable.Label);

    num_contacts = numel(targetlist);
    found_idx = nan(num_contacts, 1);

    % Match each requested target to the contact table
    for k = 1:num_contacts
        match = find(labels == targetlist(k), 1, 'first'); % use first match if duplicated
        if ~isempty(match)
            found_idx(k) = match;
        end
    end

    % Preallocate with NaNs so missing contacts stay NaN
    c2c_distance = nan(num_contacts, num_contacts);

    % Fill only pairs where both contacts were found
    for i = 1:num_contacts
        if isnan(found_idx(i))
            continue;
        end
        xi = contacttable.Coord_x(found_idx(i));
        yi = contacttable.Coord_y(found_idx(i));
        zi = contacttable.Coord_z(found_idx(i));

        for j = 1:num_contacts
            if isnan(found_idx(j))
                continue;
            end

            xj = contacttable.Coord_x(found_idx(j));
            yj = contacttable.Coord_y(found_idx(j));
            zj = contacttable.Coord_z(found_idx(j));

            c2c_distance(i,j) = sqrt((xi-xj).^2 + (yi-yj).^2 + (zi-zj).^2);
        end
    end
end