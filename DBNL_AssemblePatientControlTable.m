function [Predict, Control_Table] = DBNL_AssemblePatientControlTable(results, patientinfo, patientID)
    % Extract patient-specific data
    patientrow = find(patientinfo.sub_label == patientID);
    load(patientinfo.Montage_Location{patientrow});
    [Medlist, Medicationbinary] = DBNL_MedicationPresence(patientinfo, patientrow);
    electrodecsv = readtable(patientinfo.YAEL_Location{patientrow});
    reference_scheme = "CARLA";
    [Channel_Names, suppress_idx] = DBNL_Update_Ref(reference_scheme, Channel_Names, patientinfo, patientrow);

    Subj_Controls = [Medicationbinary, ...
                     patientinfo.age(patientrow), ...
                     patientinfo.sex{patientrow}, ...
                     patientinfo.contacts(patientrow), ...
                     patientinfo.trajectories(patientrow), ...
                     patientinfo.time_ictal(patientrow), ...
                     patientinfo.time_medication(patientrow), ...
                     patientinfo.time_reload(patientrow), ...
                     patientinfo.time_EMU(patientrow), ...
                     ['DBNLP0' num2str(patientinfo.sub_label(patientrow))], ...
                     patientinfo.Dist_ANT_to_Stim(patientrow)];

Subj_Labels = [Medlist', 'Age', 'Gender', 'num_cont', 'num_traj', 'time_ictal', 'time_med', 'time_reload', 'time_EMU', 'PatientID', 'Dist2MTT'];
                 
    % Extract dimensions of CCEP heatmap
    CCEP_heatmap = results.CCEP_Heatmap;
    CCEP_sig = ~isnan(CCEP_heatmap) & CCEP_heatmap ~= 0; % Binary significant responses
    [x, y, z] = size(CCEP_sig);
    subsample_timepoints = (0.015:0.001:0.5) * 1000;
    
    % Define time bins
    edges = [0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500];
    timeLabels = {'50ms', '100ms', '150ms', '200ms', '250ms', '300ms', '350ms', '400ms', '450ms', '500ms'};

    % Preallocate arrays
    numElements = x * y * z;
    Predict = NaN(numElements, 1);
    Control_Table_temp = cell(numElements,33); % Updated for additional column

    % Initialize valid entry counter
    valid_idx = 1;

    % Nested loops
    for i = 1:x
        if ismember(i, suppress_idx)
            % Skip suppressed rows
            continue;
        end
        parse_contact = Channel_Names(i);
        [Contact_Controls, Contact_Control_labels] = DBNL_AssembleContactControls(parse_contact, ...
            patientinfo.stim_contacts{patientrow}, ...
            patientinfo.SOZ{patientrow}, ...
            electrodecsv);

        for j = 1:y
            time_response = subsample_timepoints(j);

            % Assign time to bins
            timeBin = discretize(time_response, edges, 'categorical', timeLabels);

            for k = 1:z
                [Stim_Controls, Stim_Control_Labels] = DBNL_AssembleStimControls(patientinfo.stim_order{patientrow}, k);

                % Populate arrays with valid data
                Predict(valid_idx) = CCEP_sig(i, j, k);
                temp = [num2cell(Subj_Controls), num2cell(Contact_Controls), ...
                        {time_response, char(timeBin)}, num2cell(Stim_Controls)];
                % Flatten elements
                processedRow = cell(size(temp));
                for ii = 1:numel(temp)
                    if iscell(temp{ii}) && numel(temp{ii}) == 1
                        processedRow{ii} = temp{ii}{1};
                    else
                        processedRow{ii} = temp{ii};
                    end
                end

                % Assign the processed row
                Control_Table_temp(valid_idx, :) = processedRow;
                valid_idx = valid_idx + 1;
            end
        end
    end

    % Truncate unused rows
    Predict = Predict(1:valid_idx-1);
    Control_Table_temp = Control_Table_temp(1:valid_idx-1, :);

    % Update labels
    Control_Labels = [Subj_Labels, Contact_Control_labels, 'time', 'time_bin', Stim_Control_Labels];

    % Convert to MATLAB table
    Control_Table = cell2table(Control_Table_temp, 'VariableNames', Control_Labels);
    Control_Table.Gender = categorical(Control_Table.Gender);
    Control_Table.PatientID = categorical(Control_Table.PatientID);
    Control_Table.Cont_struct = categorical(Control_Table.Cont_struct);
    Control_Table.time_bin = categorical(Control_Table.time_bin);
    Control_Table.Pattern = categorical(Control_Table.Pattern);
    Control_Table.Current = categorical(Control_Table.Current);
    Control_Table.Cont_struct = renamecats(Control_Table.Cont_struct, 'WhiteMatter', 'Other');
    Control_Table.Sided = categorical( ...
    Control_Table.Stim_hemi == Control_Table.Cont_hemi, ...
    [true, false], ...
    {'Ipsil', 'Contra'} ...
    );

end