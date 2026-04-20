function [pairTbl, summaryTbl, out] = DBNL_deltaDistanceDependence_fromCells(resultsStruct, distMat, channelNames, varargin)
% Quantify distance-dependent similarity of pre-post difference waveforms.
% Updated to:
%   1) allow same-lead OR different-lead comparisons
%   2) compute delta-waveform peak amplitude and peak latency per contact
%   3) summarize pairwise differences in those peak features
%
% INPUTS
%   resultsStruct : patient struct with CCEP_matrix
%   distMat       : [C x C] true inter-contact distances in mm
%   channelNames  : [C x 1] labels like "RACG5"
%
% NAME-VALUE OPTIONS
%   'PatientID'      : string
%   'PairMap'        : [nBlocks x 2], PRE/POST cell mapping (default [1 2; 3 4; 5 6; 7 8])
%   'TimeVecMs'      : [1 x T] time vector in ms
%   'WindowMs'       : [start end] ms analysis window (default [15 150])
%   'TimeIdx'        : explicit indices instead of TimeVecMs/WindowMs
%   'PairMode'       : 'sameLead' (default), 'differentLead', or 'all'
%   'Metric'         : 'r2' (default), 'r', 'absr', 'maxxcorr'
%   'MaxLagSamples'  : for 'maxxcorr' (default 5)
%   'DistanceEdges'  : distance bin edges in mm
%   'SigMaskByPair'  : [C x nBlocks] logical, optional
%   'SigMode'        : 'all' (default), 'eitherSig', 'bothSig'
%   'UseAbsDelta'    : false (default)
%
% OUTPUTS
%   pairTbl    : one row per contact pair
%   summaryTbl : one row per block/pairType/distance bin
%   out        : raw outputs

%% Parse options
p = inputParser;
p.addParameter('PatientID', "", @(x)ischar(x) || isstring(x));
p.addParameter('PairMap', [1 2; 3 4; 5 6; 7 8], @(x)isnumeric(x) && size(x,2)==2);
p.addParameter('TimeVecMs', [], @(x)isnumeric(x) && isvector(x));
p.addParameter('WindowMs', [15 150], @(x)isnumeric(x) && numel(x)==2);
p.addParameter('TimeIdx', [], @(x)isnumeric(x) && isvector(x));
p.addParameter('PairMode', 'sameLead', @(x)ischar(x) || isstring(x));
p.addParameter('Metric', 'r2', @(x)ischar(x) || isstring(x));
p.addParameter('MaxLagSamples', 5, @(x)isnumeric(x) && isscalar(x));
p.addParameter('UseApproxSpacing', true, @(x)islogical(x) && isscalar(x));
p.addParameter('BaseSpacingMM', 3.5, @(x)isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('UseApproxTolerance', false, @(x)islogical(x) && isscalar(x));
p.addParameter('SpacingToleranceMM', 1.75, @(x)isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('MaxApproxSpacingMM', 70, @(x)isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('DistanceEdges', [0 5 10 20 30 40 60 inf], @(x)isnumeric(x) && isvector(x)); % fallback only
p.addParameter('SigMaskByPair', [], @(x)islogical(x) || isnumeric(x));
p.addParameter('SigMode', 'all', @(x)ischar(x) || isstring(x));
p.addParameter('UseAbsDelta', false, @(x)islogical(x) && isscalar(x));
p.parse(varargin{:});
opt = p.Results;

if isempty(opt.SpacingToleranceMM)
    opt.SpacingToleranceMM = opt.BaseSpacingMM / 2;
end

opt.PairMode = lower(string(opt.PairMode));
opt.Metric   = lower(string(opt.Metric));
opt.SigMode  = lower(string(opt.SigMode));

%% Checks
assert(isfield(resultsStruct, 'CCEP_matrix'), 'resultsStruct must contain CCEP_matrix');

C = size(distMat,1);
assert(size(distMat,2)==C, 'distMat must be C x C');
assert(numel(channelNames)==C, 'channelNames must have length C');

channelNames = string(channelNames(:));

% Parse lead and contact number from labels like RACG5
leadLabels = strings(C,1);
contactNum = nan(C,1);
for i = 1:C
    tok = regexp(channelNames(i), '^(.+?)(\d+)$', 'tokens', 'once');
    if ~isempty(tok)
        leadLabels(i) = string(tok{1});
        contactNum(i) = str2double(tok{2});
    else
        leadLabels(i) = channelNames(i);
        contactNum(i) = NaN;
        warning('Could not parse lead/contact number from channel label: %s', channelNames(i));
    end
end

nCells = numel(resultsStruct.CCEP_matrix);
assert(all(opt.PairMap(:) >= 1 & opt.PairMap(:) <= nCells), 'PairMap references invalid CCEP_matrix cells');

if isempty(opt.TimeIdx)
    assert(~isempty(opt.TimeVecMs), 'Provide TimeVecMs or TimeIdx');
    opt.TimeIdx = find(opt.TimeVecMs >= opt.WindowMs(1) & opt.TimeVecMs <= opt.WindowMs(2));
    assert(~isempty(opt.TimeIdx), 'No samples found in requested window');
end
timeWinMs = opt.TimeVecMs(opt.TimeIdx);

nBlocks = size(opt.PairMap,1);

if ~isempty(opt.SigMaskByPair)
    assert(size(opt.SigMaskByPair,1)==C, 'SigMaskByPair must be C x nBlocks');
    assert(size(opt.SigMaskByPair,2)==nBlocks, 'SigMaskByPair must match number of block pairs');
end

%% Pair masks
upperMask     = triu(true(C), 1);
validDistMask = isfinite(distMat) & ~isnan(distMat) & distMat > 0;
sameLeadMask  = leadLabels == leadLabels.';
diffLeadMask  = ~sameLeadMask;

pairMask = upperMask & validDistMask;
switch opt.PairMode
    case "samelead"
        pairMask = pairMask & sameLeadMask;
    case "differentlead"
        pairMask = pairMask & diffLeadMask;
    case "all"
        % keep all valid upper-triangle pairs
    otherwise
        error('Unknown PairMode: %s', opt.PairMode);
end

[iAll, jAll] = find(pairMask);
nBasePairs = numel(iAll);

%% Outputs
pairRows    = cell(nBlocks,1);
summaryRows = cell(nBlocks,1);

out = struct();
out.patientID = string(opt.PatientID);
out.options   = opt;
out.blockData = cell(nBlocks,1);

%% Main loop
for b = 1:nBlocks
    preIdx  = opt.PairMap(b,1);
    postIdx = opt.PairMap(b,2);

    preMat  = resultsStruct.CCEP_matrix{preIdx};   % C x T x Npre
    postMat = resultsStruct.CCEP_matrix{postIdx};  % C x T x Npost

    assert(ndims(preMat)==3 && ndims(postMat)==3, 'Each CCEP_matrix cell must be C x T x N');
    assert(size(preMat,1)==C && size(postMat,1)==C, 'First dimension must match number of contacts');

    % Average over repetitions
    preMean  = mean(preMat,  3, 'omitnan');
    postMean = mean(postMat, 3, 'omitnan');

    % Delta waveform
    deltaMean = postMean - preMean;
    if opt.UseAbsDelta
        deltaMean = abs(deltaMean);
    end

    % Restrict to analysis window
    Xw = deltaMean(:, opt.TimeIdx);

    % Per-contact peak features
    peakAmp   = nan(C,1);
    peakLatMs = nan(C,1);
    for c = 1:C
        x = Xw(c,:);
        if all(isnan(x))
            continue;
        end
        [~, idxPk] = max(abs(x));     % absolute maximum
        peakAmp(c)   = x(idxPk);      % signed amplitude at absolute max
        peakLatMs(c) = timeWinMs(idxPk);
    end

    % Optional significance filtering
    if isempty(opt.SigMaskByPair)
        sigVec = true(C,1);
    else
        sigVec = logical(opt.SigMaskByPair(:, b));
    end

    keepPair = true(nBasePairs,1);
    switch opt.SigMode
        case "all"
            % no filtering
        case "eithersig"
            keepPair = sigVec(iAll) | sigVec(jAll);
        case "bothsig"
            keepPair = sigVec(iAll) & sigVec(jAll);
        otherwise
            error('Unknown SigMode: %s', opt.SigMode);
    end

    ii = iAll(keepPair);
    jj = jAll(keepPair);

    distances = distMat(sub2ind([C C], ii, jj));
    waveformMetric = nan(numel(ii),1);

    % Assign each pair to nearest nominal electrode spacing
% Intended mainly for same-lead analyses
if opt.UseApproxSpacing
    approxDistMM = nan(numel(ii),1);
    approxStep   = nan(numel(ii),1);

    % Option A: use user-supplied base spacing directly
    baseSpacing = opt.BaseSpacingMM;

    % Option B (recommended): estimate empirical base spacing from adjacent same-lead pairs
    % Uncomment this block if you want the function to infer spacing automatically:
    %
    % if isnan(baseSpacing) || isempty(baseSpacing)
    %     adjMask = (leadLabels(ii) == leadLabels(jj)) & (abs(contactNum(ii) - contactNum(jj)) == 1);
    %     if any(adjMask)
    %         baseSpacing = median(distances(adjMask), 'omitnan');
    %     else
    %         error('Could not estimate base spacing from adjacent same-lead pairs.');
    %     end
    % end

    nominalCenters = baseSpacing:baseSpacing:opt.MaxApproxSpacingMM;

    diffToNominal = abs(distances - nominalCenters);   % implicit expansion
    [minDiff, idxNearest] = min(diffToNominal, [], 2);

    % If you do not want to impose a tolerance, assign every pair
    approxDistMM = nominalCenters(idxNearest).';
    approxStep   = idxNearest;

    % Optional QC filter if desired
    if opt.UseApproxTolerance
        keepApprox = minDiff <= opt.SpacingToleranceMM;
        approxDistMM(~keepApprox) = NaN;
        approxStep(~keepApprox)   = NaN;
    end
else
    approxDistMM = nan(numel(ii),1);
    approxStep   = nan(numel(ii),1);
end

    for k = 1:numel(ii)
        x = Xw(ii(k), :);
        y = Xw(jj(k), :);
        waveformMetric(k) = local_computeMetric(x, y, opt.Metric, opt.MaxLagSamples);
    end

    pairType = repmat("differentLead", numel(ii), 1);
    pairType(leadLabels(ii) == leadLabels(jj)) = "sameLead";

    absPeakAmpDiff = abs(peakAmp(ii) - peakAmp(jj));
    absPeakLatDiff = abs(peakLatMs(ii) - peakLatMs(jj));

    Tpair = table( ...
        repmat(string(opt.PatientID), numel(ii), 1), ...
        repmat(b, numel(ii), 1), ...
        repmat(preIdx, numel(ii), 1), ...
        repmat(postIdx, numel(ii), 1), ...
        ii, jj, ...
        channelNames(ii), channelNames(jj), ...
        leadLabels(ii), leadLabels(jj), ...
        contactNum(ii), contactNum(jj), ...
        pairType, distances, approxDistMM, approxStep, waveformMetric, ...
        peakAmp(ii), peakAmp(jj), ...
        peakLatMs(ii), peakLatMs(jj), ...
        absPeakAmpDiff, absPeakLatDiff, ...
        sigVec(ii), sigVec(jj), ...
        'VariableNames', {'PatientID','BlockPair','PreCell','PostCell', ...
        'Chan1','Chan2','ChanName1','ChanName2','Lead1','Lead2', ...
        'ContactNum1','ContactNum2','PairType','DistanceMM','ApproxDistMM','ApproxStep','WaveformMetric', ...
        'PeakAmp1','PeakAmp2','PeakLatMs1','PeakLatMs2', ...
        'AbsPeakAmpDiff','AbsPeakLatDiff','Sig1','Sig2'});

    Tpair.ContactStepSep = abs(Tpair.ContactNum1 - Tpair.ContactNum2);
    pairRows{b} = Tpair;

% Summarize by approximate nominal spacing if requested,
% otherwise fall back to broad distance edges
if opt.UseApproxSpacing
    Tvalid = Tpair(~isnan(Tpair.ApproxDistMM), :);

    if isempty(Tvalid)
        Tsum = table();
    else
        approxVals = unique(Tvalid.ApproxDistMM);
        nBins = numel(approxVals);

        medWave = nan(nBins,1);
        meanWave = nan(nBins,1);
        medAmpDiff = nan(nBins,1);
        medLatDiff = nan(nBins,1);
        nPairs = zeros(nBins,1);

        for bb = 1:nBins
            idx = Tvalid.ApproxDistMM == approxVals(bb);
            medWave(bb)    = median(Tvalid.WaveformMetric(idx), 'omitnan');
            meanWave(bb)   = mean(Tvalid.WaveformMetric(idx), 'omitnan');
            medAmpDiff(bb) = median(Tvalid.AbsPeakAmpDiff(idx), 'omitnan');
            medLatDiff(bb) = median(Tvalid.AbsPeakLatDiff(idx), 'omitnan');
            nPairs(bb)     = sum(idx);
        end

        if opt.PairMode == "all"
            pairTypeSummary = repmat("mixed", nBins, 1);
        else
            pairTypeSummary = repmat(string(opt.PairMode), nBins, 1);
        end

        Tsum = table( ...
            repmat(string(opt.PatientID), nBins, 1), ...
            repmat(b, nBins, 1), ...
            pairTypeSummary, ...
            (1:nBins)', ...
            approxVals(:), ...
            medWave, meanWave, medAmpDiff, medLatDiff, nPairs, ...
            'VariableNames', {'PatientID','BlockPair','PairType','DistanceBin','DistanceBinCenterMM', ...
            'MedianWaveformMetric','MeanWaveformMetric','MedianAbsPeakAmpDiff','MedianAbsPeakLatDiff','NPairs'});
    end

else
    Tpair.DistanceBin = discretize(Tpair.DistanceMM, opt.DistanceEdges, 'IncludedEdge', 'right');

    centers = local_edgeCenters(opt.DistanceEdges);
    nBins   = numel(centers);

    medWave = nan(nBins,1);
    meanWave = nan(nBins,1);
    medAmpDiff = nan(nBins,1);
    medLatDiff = nan(nBins,1);
    nPairs = zeros(nBins,1);

    for bb = 1:nBins
        idx = Tpair.DistanceBin == bb;
        if any(idx)
            medWave(bb)    = median(Tpair.WaveformMetric(idx), 'omitnan');
            meanWave(bb)   = mean(Tpair.WaveformMetric(idx), 'omitnan');
            medAmpDiff(bb) = median(Tpair.AbsPeakAmpDiff(idx), 'omitnan');
            medLatDiff(bb) = median(Tpair.AbsPeakLatDiff(idx), 'omitnan');
            nPairs(bb)     = sum(idx);
        end
    end

    if opt.PairMode == "all"
        pairTypeSummary = repmat("mixed", nBins, 1);
    else
        pairTypeSummary = repmat(string(opt.PairMode), nBins, 1);
    end

    Tsum = table( ...
        repmat(string(opt.PatientID), nBins, 1), ...
        repmat(b, nBins, 1), ...
        pairTypeSummary, ...
        (1:nBins)', centers(:), ...
        medWave, meanWave, medAmpDiff, medLatDiff, nPairs, ...
        'VariableNames', {'PatientID','BlockPair','PairType','DistanceBin','DistanceBinCenterMM', ...
        'MedianWaveformMetric','MeanWaveformMetric','MedianAbsPeakAmpDiff','MedianAbsPeakLatDiff','NPairs'});
end

    summaryRows{b} = Tsum;

    out.blockData{b} = struct( ...
        'BlockPair', b, ...
        'PreCell', preIdx, ...
        'PostCell', postIdx, ...
        'PreMean', preMean, ...
        'PostMean', postMean, ...
        'DeltaMean', deltaMean, ...
        'WindowedDelta', Xw, ...
        'PeakAmp', peakAmp, ...
        'PeakLatMs', peakLatMs, ...
        'PairTable', Tpair, ...
        'SummaryTable', Tsum);
end

pairTbl    = vertcat(pairRows{:});
summaryTbl = vertcat(summaryRows{:});

end


%% ========================= HELPERS ========================= %%
function m = local_computeMetric(x, y, metricName, maxLagSamples)
x = x(:);
y = y(:);

good = ~(isnan(x) | isnan(y) | isinf(x) | isinf(y));
x = x(good);
y = y(good);

if numel(x) < 5
    m = NaN;
    return;
end

switch lower(metricName)
    case 'r'
        if std(x)==0 || std(y)==0
            m = NaN;
        else
            m = corr(x, y, 'type', 'Pearson');
        end

    case 'absr'
        if std(x)==0 || std(y)==0
            m = NaN;
        else
            m = abs(corr(x, y, 'type', 'Pearson'));
        end

    case 'r2'
        if std(x)==0 || std(y)==0
            m = NaN;
        else
            r = corr(x, y, 'type', 'Pearson');
            m = r.^2;
        end

    case 'maxxcorr'
        if std(x)==0 || std(y)==0
            m = NaN;
        else
            xz = zscore(x);
            yz = zscore(y);
            c = xcorr(xz, yz, maxLagSamples, 'coeff');
            m = max(abs(c));
        end

    case 'simnrd'
        nx = norm(x, 2);
        ny = norm(y, 2);

        if nx == 0 && ny == 0
            m = NaN;
        else
            m = 1 - norm(x - y, 2) / (nx + ny);
            m = max(min(m, 1), -1);
        end

    otherwise
        error('Unknown metric: %s', metricName);
end
end

function centers = local_edgeCenters(edges)
centers = nan(numel(edges)-1,1);
for i = 1:numel(edges)-1
    if isinf(edges(i+1))
        if i == 1
            centers(i) = edges(i) + 1;
        else
            w = edges(i) - edges(i-1);
            centers(i) = edges(i) + w/2;
        end
    else
        centers(i) = mean(edges(i:i+1));
    end
end
end