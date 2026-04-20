function [heatmap_label,heatmap_idx] = DBNL_correctheatmap(sorted_FS_label)
    heatmap_label = unique(sorted_FS_label);
    for i = 1:numel(heatmap_label)
        label_count(i) = numel(find(strcmp(sorted_FS_label,heatmap_label{i})));
    end
    temp_idx = 0;
    for j = 1:numel(heatmap_label)
        heatmap_idx(j) = label_count(j)/2 + temp_idx;
        temp_idx = sum(label_count(1:j))+0.5;
    end
end