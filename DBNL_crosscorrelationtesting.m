%% Cohort local-dependence analysis:
%  A) Response window (15-300 ms), all contacts
%  B) Baseline null window (-300 to -15 ms), all contacts
%  C) Response window (15-300 ms), high-delta contacts only
%
% Metric: simNRD
% Pair type: same-lead only
% Approximate spacing bins based on true distance

%% Settings
patientIDs = ["DBNLP019","DBNLP021","DBNLP022","DBNLP023", ...
              "DBNLP025","DBNLP026","DBNLP030","DBNLP034","DBNLP035"];

PAIR_MAP         = [1 2; 3 4; 5 6; 7 8];
RESPONSE_WIN_MS  = [15 300];
BASELINE_WIN_MS  = [-300 -15];
BASE_SPACING_MM  = 3.5;
LOCAL_MAX_MM     = 21;
MAX_SPACING_MM   = 80;
HIGHDELTA_PCT    = 75;   % top 25% by delta RMS

%% Run the 3 conditions
[respPatientAcross, respCohort, respTrend, respLog] = local_collect_condition( ...
    patientIDs, RESPONSE_WIN_MS, "Response_all", false, RESPONSE_WIN_MS, ...
    PAIR_MAP, BASE_SPACING_MM, LOCAL_MAX_MM, MAX_SPACING_MM, HIGHDELTA_PCT);

[basePatientAcross, baseCohort, baseTrend, baseLog] = local_collect_condition( ...
    patientIDs, BASELINE_WIN_MS, "Baseline_all", false, RESPONSE_WIN_MS, ...
    PAIR_MAP, BASE_SPACING_MM, LOCAL_MAX_MM, MAX_SPACING_MM, HIGHDELTA_PCT);

[highPatientAcross, highCohort, highTrend, highLog] = local_collect_condition( ...
    patientIDs, RESPONSE_WIN_MS, "Response_highDelta", true, RESPONSE_WIN_MS, ...
    PAIR_MAP, BASE_SPACING_MM, LOCAL_MAX_MM, MAX_SPACING_MM, HIGHDELTA_PCT);

%% Display logs
disp('================ RESPONSE ALL LOG ================')
disp(respLog)

disp('================ BASELINE ALL LOG ================')
disp(baseLog)

disp('================ RESPONSE HIGH-DELTA LOG ================')
disp(highLog)

%% Display cohort summaries
disp('================ RESPONSE ALL COHORT SUMMARY ================')
disp(respCohort)

disp('================ BASELINE ALL COHORT SUMMARY ================')
disp(baseCohort)

disp('================ RESPONSE HIGH-DELTA COHORT SUMMARY ================')
disp(highCohort)

%% Display patient-level trend summaries
disp('================ RESPONSE ALL PATIENT TRENDS ================')
disp(respTrend)

disp('================ BASELINE ALL PATIENT TRENDS ================')
disp(baseTrend)

disp('================ RESPONSE HIGH-DELTA PATIENT TRENDS ================')
disp(highTrend)

%% Figure 1: Response vs Baseline cohort medians with patient curves
figure;
subplot(1,2,1); hold on;

% individual patient curves for response
uPatients = unique(respPatientAcross.PatientID);
for i = 1:numel(uPatients)
    T = respPatientAcross(respPatientAcross.PatientID == uPatients(i), :);
    T = sortrows(T, 'ApproxDistMM');
    plot(T.ApproxDistMM, T.PatientMedian, '-o', ...
        'Color', [0.80 0.80 0.80], 'LineWidth', 1, 'MarkerSize', 5, ...
        'HandleVisibility', 'off');
end

% cohort response
plot(respCohort.ApproxDistMM, respCohort.CohortMedian, '-o', ...
    'LineWidth', 2.5, 'MarkerSize', 7, 'DisplayName', 'Response');
for i = 1:height(respCohort)
    plot([respCohort.ApproxDistMM(i) respCohort.ApproxDistMM(i)], ...
         [respCohort.CohortQ25(i) respCohort.CohortQ75(i)], '-', ...
         'LineWidth', 1.5, 'HandleVisibility', 'off');
end

% cohort baseline
plot(baseCohort.ApproxDistMM, baseCohort.CohortMedian, '-s', ...
    'LineWidth', 2.5, 'MarkerSize', 7, 'DisplayName', 'Baseline null');
for i = 1:height(baseCohort)
    plot([baseCohort.ApproxDistMM(i) baseCohort.ApproxDistMM(i)], ...
         [baseCohort.CohortQ25(i) baseCohort.CohortQ75(i)], '-', ...
         'LineWidth', 1.5, 'HandleVisibility', 'off');
end

xlabel('Approximate inter-contact distance (mm)');
ylabel('Patient-level median simNRD');
title('Response vs baseline null');
grid on;
legend('Location', 'best');

%% Figure 2: Response all vs high-delta
subplot(1,2,2); hold on;

% cohort response all
plot(respCohort.ApproxDistMM, respCohort.CohortMedian, '-o', ...
    'LineWidth', 2.5, 'MarkerSize', 7, 'DisplayName', 'Response all');
for i = 1:height(respCohort)
    plot([respCohort.ApproxDistMM(i) respCohort.ApproxDistMM(i)], ...
         [respCohort.CohortQ25(i) respCohort.CohortQ75(i)], '-', ...
         'LineWidth', 1.5, 'HandleVisibility', 'off');
end

% cohort response high-delta
plot(highCohort.ApproxDistMM, highCohort.CohortMedian, '-^', ...
    'LineWidth', 2.5, 'MarkerSize', 7, 'DisplayName', 'Response high-delta');
for i = 1:height(highCohort)
    plot([highCohort.ApproxDistMM(i) highCohort.ApproxDistMM(i)], ...
         [highCohort.CohortQ25(i) highCohort.CohortQ75(i)], '-', ...
         'LineWidth', 1.5, 'HandleVisibility', 'off');
end

xlabel('Approximate inter-contact distance (mm)');
ylabel('Patient-level median simNRD');
title('Response all vs response high-delta');
grid on;
legend('Location', 'best');

%% Paired tests: response vs baseline at each spacing
disp('================ PAIRED TESTS: RESPONSE > BASELINE ================')
keyBins = respCohort.ApproxDistMM(:)';

for k = 1:numel(keyBins)
    d = keyBins(k);

    T1 = respPatientAcross(respPatientAcross.ApproxDistMM == d, {'PatientID','PatientMedian'});
    T2 = basePatientAcross(basePatientAcross.ApproxDistMM == d, {'PatientID','PatientMedian'});

    T1.Properties.VariableNames{'PatientMedian'} = 'RespVal';
    T2.Properties.VariableNames{'PatientMedian'} = 'BaseVal';

    J = innerjoin(T1, T2, 'Keys', 'PatientID');

    if height(J) >= 4
        p_pair = signrank(J.RespVal, J.BaseVal, 'tail', 'right');
        fprintf('Response > baseline at %.1f mm: n=%d, p=%.4g, median diff=%.4f\n', ...
            d, height(J), p_pair, median(J.RespVal - J.BaseVal, 'omitnan'));
    end
end

%% Group trend tests
validRespRho = respTrend.SpearmanRho(~isnan(respTrend.SpearmanRho));
validBaseRho = baseTrend.SpearmanRho(~isnan(baseTrend.SpearmanRho));
validHighRho = highTrend.SpearmanRho(~isnan(highTrend.SpearmanRho));

if ~isempty(validRespRho)
    p_resp = signrank(validRespRho, 0, 'tail', 'left');
    fprintf('Response-all cohort trend: median rho = %.3f, p = %.4g\n', ...
        median(validRespRho, 'omitnan'), p_resp);
end

if ~isempty(validBaseRho)
    p_base = signrank(validBaseRho, 0, 'tail', 'left');
    fprintf('Baseline-all cohort trend: median rho = %.3f, p = %.4g\n', ...
        median(validBaseRho, 'omitnan'), p_base);
end

if ~isempty(validHighRho)
    p_high = signrank(validHighRho, 0, 'tail', 'left');
    fprintf('Response-highDelta cohort trend: median rho = %.3f, p = %.4g\n', ...
        median(validHighRho, 'omitnan'), p_high);
end

%% ---------- local helper function ----------
function [patientAcross, cohortSummary, trendTbl, runLog] = local_collect_condition( ...
    patientIDs, windowMS, conditionLabel, useHighDeltaMask, highDeltaWindowMS, ...
    pairMap, baseSpacingMM, localMaxMM, maxSpacingMM, highDeltaPct)

patientAcross = table();
trendTbl      = table();
runLog        = table();

for p = 1:numel(patientIDs)

    pid = patientIDs(p);

    resultsVar = sprintf('%s_results', pid);
    distVar    = sprintf('c2c_distance_%s', pid);
    chanVar    = sprintf('%s_Channel_Names', pid);
    tVar       = sprintf('%s_t_ms', pid);

    neededVars = string({resultsVar, distVar, chanVar, tVar});

    missing = false(size(neededVars));
    for k = 1:numel(neededVars)
        missing(k) = ~evalin('base', sprintf('exist(''%s'',''var'')', neededVars(k)));
    end

    if any(missing)
        newRow = table(pid, "missing_vars", string(strjoin(neededVars(missing), ", ")), ...
            'VariableNames', {'PatientID','Status','Details'});
        runLog = [runLog; newRow];
        warning('%s skipped: missing %s', pid, strjoin(neededVars(missing), ', '));
        continue;
    end

    resultsStruct = evalin('base', resultsVar);
    distMat       = evalin('base', distVar);
    channelNames  = evalin('base', chanVar);
    t_ms          = evalin('base', tVar);

    if ~isfield(resultsStruct, 'CCEP_matrix') || isempty(resultsStruct.CCEP_matrix)
        newRow = table(pid, "bad_results", "resultsStruct.CCEP_matrix missing/empty", ...
            'VariableNames', {'PatientID','Status','Details'});
        runLog = [runLog; newRow];
        continue;
    end

    C_dist = size(distMat,1);
    C_chan = numel(channelNames);
    C_data = size(resultsStruct.CCEP_matrix{1},1);

    if size(distMat,2) ~= C_dist || C_dist ~= C_chan || C_dist ~= C_data
        details = sprintf('size mismatch: dist=%dx%d, names=%d, data=%d', ...
            size(distMat,1), size(distMat,2), C_chan, C_data);
        newRow = table(pid, "size_mismatch", string(details), ...
            'VariableNames', {'PatientID','Status','Details'});
        runLog = [runLog; newRow];
        continue;
    end

    % Optional high-delta mask computed from response window
    sigMask = [];
    sigMode = 'all';

    if useHighDeltaMask
        nBlocks = size(pairMap,1);
        C = C_data;
        sigMask = false(C, nBlocks);

        timeIdxHD = find(t_ms >= highDeltaWindowMS(1) & t_ms <= highDeltaWindowMS(2));

        for b = 1:nBlocks
            preIdx  = pairMap(b,1);
            postIdx = pairMap(b,2);

            preMean  = mean(resultsStruct.CCEP_matrix{preIdx},  3, 'omitnan');
            postMean = mean(resultsStruct.CCEP_matrix{postIdx}, 3, 'omitnan');
            deltaMean = postMean - preMean;

            X = deltaMean(:, timeIdxHD);
            deltaRMS = sqrt(mean(X.^2, 2, 'omitnan'));

            thresh = prctile(deltaRMS, highDeltaPct);
            sigMask(:,b) = deltaRMS >= thresh & ~isnan(deltaRMS);
        end

        sigMode = 'eitherSig';
    end

    % Run pairwise analysis
    [pairTbl_tmp, ~, ~] = DBNL_deltaDistanceDependence_fromCells( ...
        resultsStruct, distMat, channelNames, ...
        'PatientID', pid, ...
        'PairMap', pairMap, ...
        'TimeVecMs', t_ms, ...
        'WindowMs', windowMS, ...
        'PairMode', 'sameLead', ...
        'Metric', 'simnrd', ...
        'UseApproxSpacing', true, ...
        'BaseSpacingMM', baseSpacingMM, ...
        'UseApproxTolerance', false, ...
        'MaxApproxSpacingMM', maxSpacingMM, ...
        'SigMaskByPair', sigMask, ...
        'SigMode', sigMode);

    if isempty(pairTbl_tmp) || height(pairTbl_tmp) == 0
        newRow = table(pid, "empty_pairs", "pairTbl_tmp empty", ...
            'VariableNames', {'PatientID','Status','Details'});
        runLog = [runLog; newRow];
        continue;
    end

    pairVars = string(pairTbl_tmp.Properties.VariableNames);

    if ~ismember("ApproxDistMM", pairVars) || ~ismember("WaveformMetric", pairVars)
        newRow = table(pid, "missing_columns", "pairTbl_tmp missing ApproxDistMM or WaveformMetric", ...
            'VariableNames', {'PatientID','Status','Details'});
        runLog = [runLog; newRow];
        continue;
    end

    Tpair = pairTbl_tmp;
    Tpair = Tpair(~isnan(Tpair.ApproxDistMM), :);
    Tpair = Tpair(~isnan(Tpair.WaveformMetric), :);
    Tpair = Tpair(Tpair.ApproxDistMM <= localMaxMM, :);

    if isempty(Tpair)
        newRow = table(pid, "no_valid_local_pairs", "No valid pairs within local range", ...
            'VariableNames', {'PatientID','Status','Details'});
        runLog = [runLog; newRow];
        continue;
    end

    % Summarize by block and spacing
    [Gblock, blockVals, distVals] = findgroups(Tpair.BlockPair, Tpair.ApproxDistMM);

    Tblock = table();
    Tblock.PatientID    = repmat(pid, numel(blockVals), 1);
    Tblock.BlockPair    = blockVals;
    Tblock.ApproxDistMM = distVals;
    Tblock.NPairs       = splitapply(@numel, Tpair.WaveformMetric, Gblock);
    Tblock.BlockMedian  = splitapply(@(x) median(x,'omitnan'), Tpair.WaveformMetric, Gblock);
    Tblock.BlockQ25     = splitapply(@(x) prctile(x,25), Tpair.WaveformMetric, Gblock);
    Tblock.BlockQ75     = splitapply(@(x) prctile(x,75), Tpair.WaveformMetric, Gblock);
    Tblock.BlockMean    = splitapply(@(x) mean(x,'omitnan'), Tpair.WaveformMetric, Gblock);

    % Summarize across blocks within patient
    [Gpat, distValsPat] = findgroups(Tblock.ApproxDistMM);

    Tp = table();
    Tp.PatientID     = repmat(pid, numel(distValsPat), 1);
    Tp.Condition     = repmat(conditionLabel, numel(distValsPat), 1);
    Tp.ApproxDistMM  = distValsPat;
    Tp.NBlocks       = splitapply(@numel, Tblock.BlockMedian, Gpat);
    Tp.PatientMedian = splitapply(@(x) median(x,'omitnan'), Tblock.BlockMedian, Gpat);
    Tp.PatientQ25    = splitapply(@(x) prctile(x,25), Tblock.BlockMedian, Gpat);
    Tp.PatientQ75    = splitapply(@(x) prctile(x,75), Tblock.BlockMedian, Gpat);
    Tp.PatientMean   = splitapply(@(x) mean(x,'omitnan'), Tblock.BlockMedian, Gpat);

    patientAcross = [patientAcross; Tp];

    % Trend across spacing
    if height(Tp) >= 3
        [rho, pval] = corr(Tp.ApproxDistMM, Tp.PatientMedian, ...
            'type', 'Spearman', 'rows', 'complete');
    else
        rho = NaN;
        pval = NaN;
    end

    newTrend = table(pid, conditionLabel, rho, pval, height(Tp), ...
        'VariableNames', {'PatientID','Condition','SpearmanRho','SpearmanP','NBins'});
    trendTbl = [trendTbl; newTrend];

    newRow = table(pid, "ok", string(sprintf('ran %d spacing bins', height(Tp))), ...
        'VariableNames', {'PatientID','Status','Details'});
    runLog = [runLog; newRow];

end

% Cohort summary
if isempty(patientAcross)
    cohortSummary = table();
    return;
end

[G2, distVals2] = findgroups(patientAcross.ApproxDistMM);

cohortSummary = table();
cohortSummary.Condition    = repmat(conditionLabel, numel(distVals2), 1);
cohortSummary.ApproxDistMM = distVals2;
cohortSummary.NPatients    = splitapply(@numel, patientAcross.PatientMedian, G2);
cohortSummary.CohortMedian = splitapply(@(x) median(x,'omitnan'), patientAcross.PatientMedian, G2);
cohortSummary.CohortQ25    = splitapply(@(x) prctile(x,25), patientAcross.PatientMedian, G2);
cohortSummary.CohortQ75    = splitapply(@(x) prctile(x,75), patientAcross.PatientMedian, G2);
cohortSummary.CohortMean   = splitapply(@(x) mean(x,'omitnan'), patientAcross.PatientMedian, G2);
cohortSummary.CohortStd    = splitapply(@(x) std(x,'omitnan'), patientAcross.PatientMedian, G2);

cohortSummary = sortrows(cohortSummary, 'ApproxDistMM');

end