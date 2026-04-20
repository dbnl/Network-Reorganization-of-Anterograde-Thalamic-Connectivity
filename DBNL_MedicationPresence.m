function [Medlist, Medicationbinary] = DBNL_MedicationPresence(patientinfo, patientID)
    medicationLists = patientinfo.medication_list;
    
    % Flatten and split all medication lists to find unique medications
    allMedications = split(join(medicationLists, ', '), ', ');
    Medlist = unique(strtrim(allMedications)); % Remove duplicates and trim whitespace
    
    % Initialize binary matrix
    numMedications = length(Medlist);
    Medicationbinary = {};
    
    % Loop through each patient and populate the binary matrix
        patientspecMedications = strtrim(split(medicationLists{patientID}, ','));
        for j = 1:numMedications
            Medicationbinary{j} = ismember(Medlist{j}, patientspecMedications);
        end
    Medlist = strcat("on_", Medlist);
end