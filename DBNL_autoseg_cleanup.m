function translatedloc = DBNL_autoseg_cleanup(Channel_Names, T)
    % Define the translation mapping
    translationMap = containers.Map({'ctx-lh-pericalcarine','ctx-rh-pericalcarine','ctx-lh-lateraloccipital','ctx-rh-lateraloccipital','Right-Inf-Lat-Vent', 'Left-Inf-Lat-Vent', 'ctx-rh-supramarginal', 'ctx-lh-supramarginal', 'ctx-lh-caudalanteriorcingulate', 'ctx-lh-caudalmiddlefrontal', 'ctx-lh-cuneus', 'ctx-lh-entorhinal', 'ctx-lh-fusiform','ctx-lh-inferiorparietal', 'ctx-lh-inferiortemporal', 'ctx-lh-isthmuscingulate', 'ctx-lh-lateralorbitofrontal', 'ctx-lh-lingual', 'ctx-lh-medialorbitofrontal', 'ctx-rh-middletemporal', 'ctx-lh-parahippocampal', 'ctx-lh-paracentral', 'ctx-lh-parsopercularis', 'ctx-lh-parsorbitalis', 'ctx-lh-parstriangularis', 'ctx-lh-postcentral', 'ctx-lh-precentral', 'ctx-lh-posteriorcingulate', 'ctx-lh-precuneus', 'ctx-lh-rostralanteriorcingulate', 'ctx-lh-rostralmiddlefrontal', 'ctx-lh-superiorfrontal', 'ctx-lh-superiorparietal', 'ctx-lh-superiortemporal', 'ctx-lh-transversetemporal', 'ctx-lh-insula', 'Left-Cerebral-White-Matter', 'Left-Hippocampus', 'Left-Thalamus-Proper*', 'Left-Amygdala', 'Unknown', 'ctx-lh-middletemporal', 'ctx-rh-caudalanteriorcingulate', 'ctx-rh-caudalmiddlefrontal', 'ctx-rh-cuneus', 'ctx-rh-entorhinal', 'ctx-rh-fusiform','ctx-rh-inferiorparietal', 'ctx-rh-inferiortemporal', 'ctx-rh-isthmuscingulate', 'ctx-rh-lateralorbitofrontal', 'ctx-rh-lingual', 'ctx-rh-medialorbitofrontal', 'ctx-rh-middletemporal', 'ctx-rh-parahippocampal', 'ctx-rh-paracentral', 'ctx-rh-parsopercularis', 'ctx-rh-parsorbitalis', 'ctx-rh-parstriangularis', 'ctx-rh-postcentral', 'ctx-rh-precentral', 'ctx-rh-posteriorcingulate', 'ctx-rh-precuneus', 'ctx-rh-rostralanteriorcingulate', 'ctx-rh-rostralmiddlefrontal', 'ctx-rh-superiorfrontal', 'ctx-rh-superiorparietal', 'ctx-rh-superiortemporal', 'ctx-rh-transversetemporal', 'ctx-rh-insula', 'Right-Cerebral-White-Matter', 'Right-Hippocampus', 'Right-Thalamus-Proper*', 'Right-Amygdala', 'ctx-rh-middletemporal', 'Left-Caudate', 'Right-Caudate'}, ...
                                    {'LO-Pericalcarine','RO-Pericalcarine','LO-LateralOccipital','RO-LateralOccipital','RZ-InfLatVent', 'LZ-InfLatVent','RP-SupraMarginal', 'LP-SupraMarginal', 'LF-CaudalAnteriorCingulate', 'LF-CaudalMiddleFrontal', 'LO-Cuneus', 'LT-Entorhinal', 'LT-Fusiform', 'LP-InferiorParietal', 'LT-InferiorTemporal', 'LT-IsthmusCingulate', 'LF-LateralOrbitoFrontal', 'LO-Lingual', 'LF-MedialOrbitoFrontal', 'LT-MiddleTemporal', 'LT-Parahippocampal', 'LP-Paracentral', 'LF-ParsOpercularis', 'LF-ParsOrbitalis', 'LF-ParsTriangularis', 'LP-PostCentral', 'LF-PreCentral', 'LP-PosteriorCingulate', 'LP-Precuneus', 'LF-RostralAnteriorCingulate', 'LF-RostralMiddleFrontal', 'LF-SuperiorFrontal', 'LP-SuperiorParietal', 'LT-SuperiorTemporal', 'LT-TransverseTemporal', 'LZ-Insula', 'LZ-WhiteMatter', 'LT-Hippocampus', 'LZ-Thalamus', 'LT-Amygdala', 'ZZ-Unknown', 'LT-MiddleTemporal', 'RF-CaudalAnteriorCingulate', 'RF-CaudalMiddleFrontal', 'RO-Cuneus', 'RT-Entorhinal', 'RT-Fusiform', 'RP-InferiorParietal', 'RT-InferiorTemporal', 'RT-IsthmusCingulate', 'RF-LateralOrbitoFrontal', 'RO-Lingual', 'RF-MedialOrbitoFrontal', 'RT-MiddleTemporal', 'RT-Parahippocampal', 'RP-Paracentral', 'RF-ParsOpercularis', 'RF-ParsOrbitalis', 'RF-ParsTriangularis', 'RP-PostCentral', 'RF-PreCentral', 'RP-PosteriorCingulate', 'RP-Precuneus', 'RF-RostralAnteriorCingulate', 'RF-RostralMiddleFrontal', 'RF-SuperiorFrontal', 'RP-SuperiorParietal', 'RT-SuperiorTemporal', 'RT-TransverseTemporal', 'RZ-Insula', 'RZ-WhiteMatter', 'RT-Hippocampus', 'RZ-Thalamus', 'RT-Amygdala','RT-MiddleTemporal', 'LZ-Caudate', 'RZ-Caudate'});
    FSlabel = T.FSLabel_aparc_DKTatlas_aseg;
    % Initialize the translated array
    translatedloc = cell(size(Channel_Names));
    
    % Iterate through each string in the input array
    for i = 1:length(Channel_Names)
        idx = find(T.Label==Channel_Names(i));
        if isempty(idx)
            originalString = 'NA';
        else
            originalString = FSlabel{idx};
        end
        % Check if the string has a corresponding translation
        if isKey(translationMap, originalString)
            translatedloc{i} = translationMap(originalString);
        else
            translatedloc{i} = originalString; % Keep the original string if no translation is found
        end
    end
end