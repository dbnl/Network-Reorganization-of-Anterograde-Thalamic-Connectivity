function [Channel_Names,suppress_idx] = DBNL_Update_Ref(ref, Channel_Names, patientinfo, patientrow)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    switch ref
        case 'BP'
            channel_names_updated = strings(length(Channel_Names),1);
            channel_names_updated(end) = 'Undefined';
            
            % Loop through the elements and concatenate them
            for i = 1:length(channel_names_updated)-1
                channel_names_updated(i) = Channel_Names(i) + "-" + Channel_Names(i+1);
            end
            Channel_Names = channel_names_updated;
            [suppress_idx, ~] = DBNL_define_supress_chan(Channel_Names, patientinfo.suppress_file_location{patientrow});
            jumpidx = DBNL_detect_bipolar_change(Channel_Names);
            suppress_idx = [suppress_idx; jumpidx];
            suppress_idx = unique(suppress_idx);
        case 'LP'
            channel_names_updated = strings(length(Channel_Names),1);
            channel_names_updated(1) = 'Undefined';
            channel_names_updated(end) = 'Undefined';
            bipolar_names = strings(length(Channel_Names),1);
            for i = 1:length(channel_names_updated)-1
                bipolar_names(i) = Channel_Names(i) + "-" + Channel_Names(i+1);
            end
            % Loop through the elements and concatenate them
            for i = 2:length(channel_names_updated)-1
                channel_names_updated(i) = Channel_Names(i) + "-(" + Channel_Names(i-1) + "+" + Channel_Names(i+1) + ")";
            end
            Channel_Names = channel_names_updated;
            [suppress_idx, ~] = DBNL_define_supress_chan(Channel_Names, patientinfo.suppress_file_location{patientrow});
            jumpidx = DBNL_detect_bipolar_change(bipolar_names);
            jumpidx_start = jumpidx + 1;
            jumpidx = [jumpidx; jumpidx_start(2:end)];
            suppress_idx = [suppress_idx; 1; jumpidx];
            suppress_idx = unique(suppress_idx);
            suppress_idx(end) = [];
        otherwise
            [suppress_idx, ~] = DBNL_define_supress_chan(Channel_Names, patientinfo.suppress_file_location{patientrow});
    end
end