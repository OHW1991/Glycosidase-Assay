function plotCoGDistanceSummaries_LAMP1(analysisRoot)
%PLOTCOGDISTANCESUMMARIES_LAMP1  Distance summaries for LAMP1 CoGs.
%
%   plotCoGDistanceSummaries_LAMP1('AnalysisDistance')
%
% Uses Excel files from:
%   AnalysisDistance/<SampleName>/*_distances.xlsx
%
% Figure 1 (2x1):
%   Row 1 – EndoH pairs, SNA- vs SNA+ (blue vs green)
%   Row 2 – NA pairs,   SNA- vs SNA+ (blue vs green)
%   For each pair:
%       - boxchart + orange swarm + black X mean
%       - annotations: p-value, ΔD = mean(SNA+) − mean(SNA−), N-/N+
%
% Figure 2 (2x2):
%   [1,1] EndoH SNA-
%   [1,2] EndoH SNA+
%   [2,1] NA    SNA-
%   [2,2] NA    SNA+
%   In each panel:
%       - one box per condition
%       - N above each box
%       - pairwise t-tests (log10 distances) between all unordered pairs,
%         shown only above the LEFT box of each pair.
%
% Figure 3:
%   ΔD trendlines vs enzyme level, using EXACTLY the ΔD values from
%   Figure 1 (so the numbers match 1:1).
%       • EndoH series @ 37°C:   EndoH−, 1:20, 1:2
%       • NA series:             NB, N10, N2
%   X-ticks:  No enzyme, Low enzyme, High enzyme.
%
% All y-axes are log10, from 1e-3 to 1e1 µm.

    if nargin < 1 || isempty(analysisRoot)
        analysisRoot = 'AnalysisDistance';
    end
    if ~isfolder(analysisRoot)
        error('AnalysisDistance root "%s" not found.', analysisRoot);
    end

    pixToUm = 101.61 / 512;  % µm per pixel

    % ---------------------------------------------------------------------
    % 1. Define the pairings (SNA- vs SNA+)
    % ---------------------------------------------------------------------

    % EndoH series – order on the x-axis:
    %   1) 37C EndoH- (R8L+EB)
    %   2) 37C EndoH 1:20 (E10)
    %   3) 37C EndoH 1:2  (R8L+E2)
    %   4) 25C EndoH 1:2  (R8L)
    endoPairs(1).label   = '37C EndoH- (R8L+EB)';
    endoPairs(1).minus   = {'R8L 37C EndoH- SNA-','EB'};
    endoPairs(1).plus    = {'R8L 37C EndoH- SNA+','EB SNA'};

    endoPairs(2).label   = '37C EndoH 1:20 (E10)';
    endoPairs(2).minus   = {'E10'};
    endoPairs(2).plus    = {'E10 SNA'};

    endoPairs(3).label   = '37C EndoH 1:2 (R8L+E2)';
    endoPairs(3).minus   = {'R8L 37C EndoH+ SNA-','E2'};
    endoPairs(3).plus    = {'R8L 37C EndoH+ SNA+','E2 SNA'};

    endoPairs(4).label   = '25C EndoH 1:2 (R8L)';
    endoPairs(4).minus   = {'R8L 25C EndoH+ SNA-'};
    endoPairs(4).plus    = {'R8L 25C EndoH+ SNA+'};

    % NA series
    naPairs(1).label     = 'NB';
    naPairs(1).minus     = {'NB'};
    naPairs(1).plus      = {'NB SNA'};

    naPairs(2).label     = 'N10';
    naPairs(2).minus     = {'N10'};
    naPairs(2).plus      = {'N10 SNA'};

    naPairs(3).label     = 'N2';
    naPairs(3).minus     = {'N2'};
    naPairs(3).plus      = {'N2 SNA'};

    % ---------------------------------------------------------------------
    % 2. Read / aggregate all distances once
    % ---------------------------------------------------------------------
    numEndo = numel(endoPairs);
    numNA   = numel(naPairs);

    dEndoMinus = cell(1,numEndo);
    dEndoPlus  = cell(1,numEndo);
    for k = 1:numEndo
        dEndoMinus{k} = readSampleListDistances(analysisRoot, endoPairs(k).minus, pixToUm);
        dEndoPlus{k}  = readSampleListDistances(analysisRoot, endoPairs(k).plus,  pixToUm);
    end

    dNAMinus = cell(1,numNA);
    dNAPlus  = cell(1,numNA);
    for k = 1:numNA
        dNAMinus{k} = readSampleListDistances(analysisRoot, naPairs(k).minus, pixToUm);
        dNAPlus{k}  = readSampleListDistances(analysisRoot, naPairs(k).plus,  pixToUm);
    end

    endoLabels = {endoPairs.label};
    naLabels   = {naPairs.label};

    % ---------------------------------------------------------------------
    % 3. Figure 1 – paired SNA- vs SNA+
    % ---------------------------------------------------------------------
    fig1 = figure('Color','w','Name','LAMP1 CoG distances: SNA- vs SNA+');
    set(fig1,'Units','normalized','OuterPosition',[0 0 1 1]);

    subplot(2,1,1); hold on;
    title('EndoH series: SNA^- vs SNA^+ (RNA^+)');
    [deltaEndoAll, ~, ~, pEndoAll] = makePairRow(dEndoMinus, dEndoPlus, endoLabels);
    ylabel('Distance (\mum)');

    subplot(2,1,2); hold on;
    title('NA series: SNA^- vs SNA^+ (RNA^+)');
    [deltaNAAll,  ~, ~, pNAAll] = makePairRow(dNAMinus,  dNAPlus,  naLabels);
    ylabel('Distance (\mum)');

    % ---------------------------------------------------------------------
    % 4. Figure 2 – grouped by SNA status (within-SNA pairwise tests)
    % ---------------------------------------------------------------------
    fig2 = figure('Color','w','Name','LAMP1 CoG distances: grouped by SNA');
    set(fig2,'Units','normalized','OuterPosition',[0 0 1 1]);

    subplot(2,2,1); hold on;
    plotGroupPanel(dEndoMinus, endoLabels, 'EndoH SNA-', [0.30 0.50 0.90]);

    subplot(2,2,2); hold on;
    plotGroupPanel(dEndoPlus,  endoLabels, 'EndoH SNA+', [0.40 0.80 0.40]);

    subplot(2,2,3); hold on;
    plotGroupPanel(dNAMinus,   naLabels,   'NA SNA-',   [0.30 0.50 0.90]);

    subplot(2,2,4); hold on;
    plotGroupPanel(dNAPlus,    naLabels,   'NA SNA+',   [0.40 0.80 0.40]);

    % ---------------------------------------------------------------------
    % 5. Figure 3 – ΔD trendlines (using EXACT ΔD & p-values from Figure 1)
    % ---------------------------------------------------------------------
    % EndoH trend: only 37°C points → indices 1,2,3
    deltaEndoTrend = deltaEndoAll(1:3);
    pEndoTrend     = pEndoAll(1:3);

    % NA trend: NB, N10, N2 → indices 1,2,3
    deltaNATrend   = deltaNAAll(1:3);
    pNATrend       = pNAAll(1:3);

    x       = 1:3;
    xLabels = {'No enzyme','Low enzyme','High enzyme'};

    % --- convert p-values to -log10(p), then normalise 0–1 ---
    clipP       = @(p) max(min(p,1), 1e-10);          % avoid 0 / >1
    logP_Endo   = -log10(clipP(pEndoTrend));
    logP_NA     = -log10(clipP(pNATrend));

    allLogP     = [logP_Endo(:); logP_NA(:)];
    allLogP     = allLogP(~isnan(allLogP));
    if isempty(allLogP)
        maxLogP = 1;
    else
        maxLogP = max(allLogP);
    end

    normEndo = logP_Endo / maxLogP;   % 0 → lowest sig, 1 → highest sig
    normNA   = logP_NA   / maxLogP;

    % sizes: area 60–200
    minSz = 60;
    maxSz = 200;
    sizeEndo = minSz + normEndo * (maxSz - minSz);
    sizeNA   = minSz + normNA   * (maxSz - minSz);

    % colours: blue → white → red through colourFromNorm()
    colorEndo = zeros(numel(x),3);
    colorNA   = zeros(numel(x),3);
    for i = 1:numel(x)
        colorEndo(i,:) = colorFromNorm(normEndo(i));
        colorNA(i,:)   = colorFromNorm(normNA(i));
    end

    % y-limits: 0 – max(ΔD)
    allDeltaTrend = [deltaEndoTrend(:); deltaNATrend(:)];
    allDeltaTrend = allDeltaTrend(~isnan(allDeltaTrend));
    if isempty(allDeltaTrend)
        yMaxTrend = 1;
    else
        yMaxTrend = max(allDeltaTrend) * 1.1;
    end

    fig3 = figure('Color','w','Name','\DeltaD trend: EndoH vs NA');
    set(fig3,'Units','normalized','OuterPosition',[0.15 0.15 0.5 0.5]);
    hold on;

    % black lines (EndoH solid, NA dashed)
    plot(x, deltaEndoTrend, 'k-',  'LineWidth', 1.5);
    plot(x, deltaNATrend,   'k--', 'LineWidth', 1.5);

    % circles: colour & size from -log10(p); slight x-offset so pairs are visible
    for i = 1:numel(x)
        if ~isnan(deltaEndoTrend(i))
            scatter(x(i)-0.03, deltaEndoTrend(i), sizeEndo(i), ...
                colorEndo(i,:), 'filled', 'MarkerEdgeColor','k');
        end
        if ~isnan(deltaNATrend(i))
            scatter(x(i)+0.03, deltaNATrend(i), sizeNA(i), ...
                colorNA(i,:), 'filled', 'MarkerEdgeColor','k');
        end
    end

    % ---------------------------------------------------------------------
    % Add colorbar for –log10(p)
    % ---------------------------------------------------------------------
    % Construct a fake colormap vector from 0→1 (the norm range)
    cmapVals = linspace(0,1,256)';
    cmapRGB  = zeros(256,3);
    for ii = 1:256
        cmapRGB(ii,:) = colorFromNorm(cmapVals(ii));
    end
    colormap(cmapRGB);

    % Determine the actual p-value range used
    allLogPTrend = [logP_Endo(:); logP_NA(:)];
    allLogPTrend = allLogPTrend(~isnan(allLogPTrend));
    if isempty(allLogPTrend)
        maxLogPTrend = 1;
    else
        maxLogPTrend = max(allLogPTrend);
    end

    cb = colorbar;
    cb.Label.String = '-log_{10}(p)';
    cb.Ticks = linspace(0,1,5);
    cb.TickLabels = arrayfun(@(v) sprintf('%.2g', v*maxLogPTrend), cb.Ticks, 'UniformOutput', false);

    % Match caxis to the normalized –log10(p)
    caxis([0 1]);

    xlim([0.8 3.2]);
    ylim([0 yMaxTrend]);
    xticks(x);
    xticklabels(xLabels);

    ylabel('\DeltaD (mean SNA^+ - mean SNA^-) [\mum]');
    grid on;
    legend({'EndoH series','NA series'}, 'Location','best');
    title('\DeltaD vs enzyme level');

    hold off;
    box on;
end

% =====================================================================
% Helper: read all distances for a *list* of sample names and pool them
% =====================================================================
function dUm = readSampleListDistances(rootDir, sampleNames, pixToUm)
    dUm = [];

    allDirs  = dir(rootDir);
    allDirs  = allDirs([allDirs.isdir]);
    allNames = {allDirs.name};

    for i = 1:numel(sampleNames)
        pattern  = sampleNames{i};
        exactDir = fullfile(rootDir, pattern);
        matchedDirs = {};

        if isfolder(exactDir)
            matchedDirs = {pattern};
        else
            patUpper = upper(pattern);
            for j = 1:numel(allNames)
                if any(strcmp(allNames{j},{'.','..'})), continue; end
                if contains(upper(allNames{j}), patUpper)
                    matchedDirs{end+1} = allNames{j}; %#ok<AGROW>
                end
            end
        end

        if isempty(matchedDirs)
            fprintf('Warning: no folder matched pattern "%s" in %s\n', ...
                pattern, rootDir);
            continue;
        end

        for j = 1:numel(matchedDirs)
            sdir = fullfile(rootDir, matchedDirs{j});
            dUm  = [dUm; readAllDistances(sdir, pixToUm)]; %#ok<AGROW>
        end
    end
end

% ---------------------------------------------------------------------
function dUm = readAllDistances(sampleDir, pixToUm)
    dUm = [];
    if ~isfolder(sampleDir), return; end
    files = dir(fullfile(sampleDir,'*_distances.xlsx'));
    for iF = 1:numel(files)
        fpath = fullfile(sampleDir, files(iF).name);
        try
            T = readtable(fpath);
        catch
            continue;
        end
        if ~ismember('DistancePixels', T.Properties.VariableNames)
            continue;
        end
        d = T.DistancePixels(:);
        d = d(~isnan(d) & d > 0);
        dUm = [dUm; d * pixToUm]; %#ok<AGROW>
    end
end

% =====================================================================
% Figure 1: paired SNA- vs SNA+ row
% =====================================================================
% =====================================================================
% Figure 1: paired SNA- vs SNA+ row
% =====================================================================
function [deltaVec, meanMinusVec, meanPlusVec, pVec] = makePairRow(dMinusCell, dPlusCell, labels)
    numPairs = numel(labels);

    % output arrays – used later for trendline & colormap
    deltaVec     = nan(1,numPairs);
    meanMinusVec = nan(1,numPairs);
    meanPlusVec  = nan(1,numPairs);
    pVec         = nan(1,numPairs);   % <- NEW

    % colours
    colMinus = [0.30 0.50 0.90];  % blue
    colPlus  = [0.40 0.80 0.40];  % green
    colSwarm = [1.00 0.65 0.20];  % orange

    logY1 = -3;  logY2 = 1;
    yP = 10^(logY1 + 0.92*(logY2-logY1));  % p-value
    yD = 10^(logY1 + 0.82*(logY2-logY1));  % ΔD
    yN = 10^(logY1 + 0.72*(logY2-logY1));  % N-/N+

    for k = 1:numPairs
        minusD = dMinusCell{k};
        plusD  = dPlusCell{k};

        xMinus = k - 0.15;
        xPlus  = k + 0.15;

        % --- means (can be NaN if empty) ---
        muMinus = mean(minusD, 'omitnan');
        muPlus  = mean(plusD,  'omitnan');
        deltaD  = muPlus - muMinus;

        meanMinusVec(k) = muMinus;
        meanPlusVec(k)  = muPlus;
        deltaVec(k)     = deltaD;

        % --- SNA- box + swarm + mean ---
        if ~isempty(minusD)
            boxchart(ones(size(minusD))*xMinus, minusD, ...
                'BoxWidth',0.25, 'BoxFaceColor',colMinus, ...
                'MarkerStyle','none', 'BoxFaceAlpha',0.5);
            swarmchart(ones(size(minusD))*xMinus, minusD, 10, ...
                'filled','MarkerFaceColor',colSwarm, ...
                'MarkerEdgeColor','none','MarkerFaceAlpha',0.6, ...
                'XJitter','density','XJitterWidth',0.15);
            plot(xMinus, muMinus, 'kx', 'MarkerSize',8, 'LineWidth',1.2);
        end

        % --- SNA+ box + swarm + mean ---
        if ~isempty(plusD)
            boxchart(ones(size(plusD))*xPlus, plusD, ...
                'BoxWidth',0.25, 'BoxFaceColor',colPlus, ...
                'MarkerStyle','none', 'BoxFaceAlpha',0.5);
            swarmchart(ones(size(plusD))*xPlus, plusD, 10, ...
                'filled','MarkerFaceColor',colSwarm, ...
                'MarkerEdgeColor','none','MarkerFaceAlpha',0.6, ...
                'XJitter','density','XJitterWidth',0.15);
            plot(xPlus, muPlus, 'kx', 'MarkerSize',8, 'LineWidth',1.2);
        end

        % --- stats & annotations ---
        if ~isempty(minusD) && ~isempty(plusD) && numel(minusD) > 1 && numel(plusD) > 1
            [~, p] = ttest2((minusD), (plusD));
            pVec(k) = p;   % <- NEW: store p-value

            text(k, yP, sprintf('p=%.2g', p), ...
                'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
                'FontSize',8);
            text(k, yD, sprintf('\\DeltaD=%.2f', deltaD), ...
                'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
                'FontSize',8);
        end

        text(k, yN, sprintf('N-/N+ = %d/%d', numel(minusD), numel(plusD)), ...
            'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
            'FontSize',8);
    end

    set(gca,'YScale','log','YLim',[1e-3 1e1]);
    xlim([0.5 numPairs+0.5]);
    xticks(1:numPairs);
    xticklabels(labels);
    xtickangle(45);
    grid on; box on;
end

% =====================================================================
% Figure 2: grouped panels for each SNA status
% =====================================================================
function plotGroupPanel(dCell, labels, panelTitle, boxColor)
    ax = gca;
    cla(ax); hold(ax,'on');

    colSwarm = [1.00 0.65 0.20];  % orange swarm
    B        = numel(labels);

    % clear appdata from previous uses
    for k = 1:B
        tagV = sprintf('vals_box_%d',k);
        tagN = sprintf('N_box_%d',k);
        if isappdata(ax,tagV), rmappdata(ax,tagV); end
        if isappdata(ax,tagN), rmappdata(ax,tagN); end
    end

    for k = 1:B
        vals = dCell{k};
        if isempty(vals), continue; end

        x = k * ones(size(vals));

        boxchart(x, vals, ...
            'BoxWidth',0.35, 'BoxFaceColor',boxColor, ...
            'MarkerStyle','none','BoxFaceAlpha',0.5);
        swarmchart(x, vals, 10, ...
            'filled','MarkerFaceColor',colSwarm, ...
            'MarkerEdgeColor','none','MarkerFaceAlpha',0.6, ...
            'XJitter','density','XJitterWidth',0.20);
        plot(k, mean(vals,'omitnan'), 'kx', 'MarkerSize',8, 'LineWidth',1.2);

        setappdata(ax, sprintf('vals_box_%d',k), vals);
        setappdata(ax, sprintf('N_box_%d',k),   numel(vals));
    end

    title(panelTitle);
    ylabel('Distance (\mum)');
    formatLogYWithN_andPairwise(labels);
end

% =====================================================================
function formatLogYWithN_andPairwise(labels)
    ax = gca;
    B  = numel(labels);

    set(ax,'YScale','log','YLim',[1e-3 1e1]);

    logY1 = -3;  logY2 = 1;
    yN        = 10^(logY1 + 0.70*(logY2-logY1));   % N line
    yTop      = 10^(logY1 + 0.92*(logY2-logY1));   % first pairwise line
    yStepLog  = 0.12*(logY2-logY1);                % spacing in log10

    % N annotations
    for k = 1:B
        tagN = sprintf('N_box_%d',k);
        if isappdata(ax,tagN)
            Nk = getappdata(ax,tagN);
            if ~isempty(Nk) && Nk > 0
                text(k, yN, sprintf('N=%d', Nk), ...
                    'HorizontalAlignment','center', ...
                    'VerticalAlignment','bottom', ...
                    'FontSize',8);
            end
        end
    end

    % pairwise tests
    if B >= 2
        for k = 1:B
            tagV_k = sprintf('vals_box_%d',k);
            if ~isappdata(ax,tagV_k), continue; end
            valsK = getappdata(ax,tagV_k);
            if numel(valsK) < 2, continue; end

            lineIdx = 0;
            for m = k+1:B
                tagV_m = sprintf('vals_box_%d',m);
                if ~isappdata(ax,tagV_m), continue; end
                valsM = getappdata(ax,tagV_m);
                if numel(valsM) < 2, continue; end

                lineIdx = lineIdx + 1;

                [~, p] = ttest2((valsK), (valsM));
                deltaD = mean(valsM,'omitnan') - mean(valsK,'omitnan');

                logY = log10(yTop) - (lineIdx-1)*yStepLog;
                yPos = 10^(logY);

                otherLabel = labels{m};
                if numel(otherLabel) > 20
                    otherLabel = [otherLabel(1:20) '…'];
                end

                txt = sprintf('%s: p=%.2g, \\DeltaD=%.2f', otherLabel, p, deltaD);
                text(k, yPos, txt, ...
                    'HorizontalAlignment','center', ...
                    'VerticalAlignment','bottom', ...
                    'FontSize',7);
            end
        end
    end

    xlim([0.5 B+0.5]);
    xticks(1:B);
    xticklabels(labels);
    xtickangle(45);
    grid on; box on;
end

% =====================================================================
% Helper: map a value in [0,1] to a blue–white–red colour
% =====================================================================
function col = colorFromNorm(t)
    if isnan(t)
        col = [0.5 0.5 0.5];  % grey for missing
        return;
    end
    t = max(0,min(1,t));      % clamp

    blue  = [0 0 1];
    white = [1 1 1];
    red   = [1 0 0];

    if t <= 0.5
        a   = t / 0.5;               % 0→blue, 0.5→white
        col = (1-a)*blue + a*white;
    else
        a   = (t-0.5) / 0.5;         % 0.5→white, 1→red
        col = (1-a)*white + a*red;
    end
end
