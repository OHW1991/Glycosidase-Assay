function plotRGandRB_counts_R8L_SNAplus(eventsRoot)
%PLOTRGANDRB_COUNTS_R8L_SNAPLUS
%   2×1 figure of RG (incl triple) and RB-only counts, only:
%       • SNA+ samples
%       • non–R-L masters (R8L + 20251130 samples)
%
%   Top row – EndoH series (SNA+, NEW ORDER):
%       1) 37C EndoH- SNA+          (R8L 37C EndoH- SNA+)
%       2) 37C EndoH 1:20 (E10 SNA) (E10 SNA)
%       3) 37C EndoH 1:2 SNA+       (R8L 37C EndoH+ SNA+)
%       4) 25C EndoH 1:2 SNA+       (R8L 25C EndoH+ SNA+)
%
%   Bottom row – NA concentration series (SNA+):
%       NB SNA
%       N10 SNA
%       N2 SNA
%
%   RG = events with ColocGreen==1 OR ColocAll==1
%   RB-only = events with ColocBlue==1 AND ColocAll==0

    if nargin < 1 || isempty(eventsRoot)
        eventsRoot = 'EventsAnalysis';
    end

    if ~isfolder(eventsRoot)
        error('Events root "%s" not found.', eventsRoot);
    end

    % Collect all *_events.xlsx and their condition keys + meta
    [eventFiles, condKeys, meta] = collectEventFiles(eventsRoot);

    % Filter to SNA+ and NOT R-L
    isSNA = contains({meta.sample}', 'SNA', 'IgnoreCase', true);
    isRL  = strcmp({meta.master}', 'R-L');
    keep  = isSNA & ~isRL;

    eventFiles = eventFiles(keep);
    condKeys   = condKeys(keep);
    meta       = meta(keep); %#ok<NASGU>

    if isempty(eventFiles)
        warning('No SNA+ non–R-L event files found for plotting.');
        return;
    end

    % Unique condition keys (after pooling via conditionKey)
    keys = unique(condKeys, 'stable');
    nK   = numel(keys);

    % Accumulate RG and RB-only counts across FoVs per conditionKey
    RG = zeros(1, nK);
    RB = zeros(1, nK);

    for i = 1:numel(eventFiles)
        key = condKeys{i};
        k   = find(strcmp(keys, key), 1);

        T = readtable(eventFiles{i});
        if isempty(T), continue; end

        CG = logical(T.ColocGreen);
        CB = logical(T.ColocBlue);
        CA = logical(T.ColocAll);

        RB(k) = RB(k) + nnz(CB & ~CA);   % RB-only
        RG(k) = RG(k) + nnz(CG | CA);    % RG + triple
    end

    % Overall max for y-limits
    allCounts  = [RG, RB];
    globalMax  = max(allCounts);
    if isempty(globalMax) || globalMax <= 0
        globalMax = 1;
    end
    globalYmax = globalMax * 1.2;

    % Common colours
    edgeCol   = [0.4 0 0];  % dark red edge
    faceGreen = [0 0.7 0];  % green fill (RG)
    faceBlue  = [0 0.4 1];  % blue fill (RB)

    % ---------- Define rows & desired orders ----------

    % Top row: EndoH series (new order: EndoH-, 1:20, 1:2 37C, 1:2 25C)
    rowEndo_labels = { ...
        '37C EndoH- SNA+', ...
        '37C EndoH 1:20 (E10 SNA)', ...
        '37C EndoH 1:2 SNA+', ...
        '25C EndoH 1:2 SNA+'};

    rowEndo_patterns = { ...
        {{'R8L','37','EndoH-','SNA+'}, {}}, ...
        {{'E10','SNA'}, {}}, ...
        {{'R8L','37','EndoH+','SNA+'}, {}}, ...
        {{'R8L','25','EndoH+','SNA+'}, {}}};

    % Bottom row: NA concentration series
    rowNA_labels = { ...
        'NB SNA', ...
        'N10 SNA', ...
        'N2 SNA'};

    rowNA_patterns = { ...
        {{'NB','SNA'}, {}}, ...
        {{'N10','SNA'}, {}}, ...
        {{'N2','SNA'}, {}}};

    % ---------- Build figure (2×1) ----------
    fig = figure('Color','w', ...
                 'Name', 'RG & RB counts (R8L+NEW, SNA+ only) – 2×1');
    tl = tiledlayout(fig, 2, 1, 'TileSpacing','compact', 'Padding','compact');

    % Top row: EndoH
    ax1 = nexttile(tl, 1);
    plotRow(ax1, rowEndo_labels, rowEndo_patterns, keys, RG, RB, ...
            faceGreen, faceBlue, edgeCol, globalYmax);
    title(ax1, 'EndoH series (R8L + E10, SNA+)');
    box(ax1, 'on');

    % Bottom row: NA
    ax2 = nexttile(tl, 2);
    plotRow(ax2, rowNA_labels, rowNA_patterns, keys, RG, RB, ...
            faceGreen, faceBlue, edgeCol, globalYmax);
    title(ax2, 'NA concentration series (NB/N10/N2, SNA+)');
    box(ax2, 'on');

    % ---------- Legend on top axes ----------
    hold(ax1, 'on');
    hRG = plot(ax1, nan, nan, 'o', ...
        'MarkerFaceColor',  faceGreen, ...
        'MarkerEdgeColor',  edgeCol, ...
        'MarkerSize',       8, ...
        'LineWidth',        1.5);
    hRB = plot(ax1, nan, nan, 'o', ...
        'MarkerFaceColor',  faceBlue, ...
        'MarkerEdgeColor',  edgeCol, ...
        'MarkerSize',       8, ...
        'LineWidth',        1.5);
    hold(ax1, 'off');

    legend(ax1, [hRG, hRB], {'RG (incl RGB)','RB-only'}, ...
        'Location','northoutside', 'Orientation','horizontal');
end

% =====================================================================
function plotRow(ax, row_labels, row_patterns, keys, RG, RB, ...
                 faceGreen, faceBlue, edgeCol, globalYmax)
% Plot a single row (one axes) of RG/RB counts

    n = numel(row_labels);
    x = 1:n;

    axes(ax); %#ok<LAXES>
    hold(ax, 'on');

    for j = 1:n
        incl = row_patterns{j}{1};
        excl = row_patterns{j}{2};
        idx  = findKeyIdxByPattern(keys, incl, excl);

        if isnan(idx)
            yRG = 0;
            yRB = 0;
        else
            yRG = RG(idx);
            yRB = RB(idx);
        end

        % RG point
        if yRG > 0
            plot(ax, x(j), yRG, 'o', ...
                'MarkerFaceColor',  faceGreen, ...
                'MarkerEdgeColor',  edgeCol, ...
                'MarkerSize',       8, ...
                'LineWidth',        1.5);
        end

        % RB point
        if yRB > 0
            plot(ax, x(j), yRB, 'o', ...
                'MarkerFaceColor',  faceBlue, ...
                'MarkerEdgeColor',  edgeCol, ...
                'MarkerSize',       8, ...
                'LineWidth',        1.5);
        end
    end

    hold(ax, 'off');

    ylim(ax, [0, 100]);%globalYmax]);
    xlim(ax, [0.5, n + 0.5]);
    xticks(ax, 1:n);
    xticklabels(ax, row_labels);
    xtickangle(ax, 45);
    ylabel(ax, 'Number of events');
end

% =====================================================================
function idx = findKeyIdxByPattern(keys, inclSubstrings, exclSubstrings)
%FINDKEYIDXBYPATTERN  Find the index in KEYS that matches all INCL and none of EXCL.
%
%   idx = findKeyIdxByPattern(keys, {'R8L','37','EndoH+','SNA+'}, {'E10'})
%
% Returns NaN if no match.

    if nargin < 3 || isempty(exclSubstrings)
        exclSubstrings = {};
    end

    mask = true(size(keys));

    % required substrings
    for i = 1:numel(inclSubstrings)
        s = inclSubstrings{i};
        mask = mask & contains(keys, s, 'IgnoreCase', true);
    end

    % forbidden substrings
    for i = 1:numel(exclSubstrings)
        s = exclSubstrings{i};
        mask = mask & ~contains(keys, s, 'IgnoreCase', true);
    end

    idxCandidates = find(mask);
    if isempty(idxCandidates)
        idx = NaN;
    else
        idx = idxCandidates(1);
    end
end
