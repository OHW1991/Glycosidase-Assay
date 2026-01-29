function plotRGoverRB_R8L_SNAplus(eventsRoot)
%PLOTRGOVERRB_R8L_SNAPLUS
%   2×1 bar-plot of RG(incl. RGB) / RB-only ratio, only:
%       • SNA+ samples
%       • non–R-L masters
%
%   Ratio is computed per conditionKey:
%       RG = ColocGreen == 1 OR ColocAll == 1
%       RB-only = ColocBlue == 1 AND ColocAll == 0
%
%   Top – EndoH series (new order):
%       1) 37C EndoH- SNA+          (R8L 37C EndoH- SNA+)
%       2) 37C EndoH 1:20 (E10 SNA) (E10 SNA)
%       3) 37C EndoH 1:2 SNA+       (R8L 37C EndoH+ SNA+)
%       4) 25C EndoH 1:2 SNA+       (R8L 25C EndoH+ SNA+)
%
%   Bottom – NA concentration series:
%       NB SNA, N10 SNA, N2 SNA

    if nargin < 1 || isempty(eventsRoot)
        eventsRoot = 'EventsAnalysis';
    end
    if ~isfolder(eventsRoot)
        error('Events root "%s" not found.', eventsRoot);
    end

    % ---- collect event files ----
    [eventFiles, condKeys, meta] = collectEventFiles(eventsRoot);

    % SNA+ only, non–R-L
    isSNA = contains({meta.sample}', 'SNA', 'IgnoreCase', true);
    isRL  = strcmp({meta.master}', 'R-L');
    keep  = isSNA & ~isRL;

    eventFiles = eventFiles(keep);
    condKeys   = condKeys(keep);
    meta       = meta(keep); %#ok<NASGU>

    if isempty(eventFiles)
        warning('No SNA+ non–R-L event files found for RG/RB ratio plot.');
        return;
    end

    % Unique condition keys (after pooling via conditionKey)
    keys = unique(condKeys, 'stable');
    nK   = numel(keys);

    % Accumulate RG and RB-only counts across FoVs
    RG = zeros(1,nK);
    RB = zeros(1,nK);

    for i = 1:numel(eventFiles)
        k = find(strcmp(keys, condKeys{i}), 1);
        T = readtable(eventFiles{i});
        if isempty(T), continue; end

        CG = logical(T.ColocGreen);
        CB = logical(T.ColocBlue);
        CA = logical(T.ColocAll);

        RB(k) = RB(k) + nnz(CB & ~CA);   % RB-only
        RG(k) = RG(k) + nnz(CG | CA);    % RG + triple
    end

    % ---- compute ratios per key ----
    ratioTrue = nan(1,nK);
    hasRB     = RB > 0;
    ratioTrue(hasRB) = RG(hasRB) ./ RB(hasRB);

    ratioPlot = zeros(1,nK);
    ratioPlot(hasRB) = ratioTrue(hasRB);

    finiteVals = ratioTrue(~isnan(ratioTrue) & isfinite(ratioTrue));
    if isempty(finiteVals)
        base = 1;
    else
        base = max(finiteVals);
        if base <= 0, base = 1;
    end

    % “Infinite” ratios: RB=0 & RG>0
    infIdx = (RB == 0) & (RG > 0);
    ratioPlot(infIdx) = base * 1.05;

    % Global y-limit
    allPos = ratioPlot(ratioPlot > 0);
    if isempty(allPos)
        ymax = 1;
    else
        ymax = max(allPos) * 1.2;
    end

    % ---------- DEFINE ROWS (2×1) ----------

    % Top row: EndoH series – order: EndoH-, 1:20, 1:2 (37C), 1:2 (25C)
    endo_labels = { ...
        '37C EndoH- SNA+', ...
        '37C EndoH 1:20 (E10 SNA)', ...
        '37C EndoH 1:2 SNA+', ...
        '25C EndoH 1:2 SNA+'};

    endo_patterns = { ...
        {{'R8L','37','EndoH-','SNA+'}, {}}, ...
        {{'E10','SNA'}, {}}, ...
        {{'R8L','37','EndoH+','SNA+'}, {}}, ...
        {{'R8L','25','EndoH+','SNA+'}, {}}};

    % Bottom row: NA series (NB, N10, N2) – like CoG figure
    na_labels = { ...
        'NB SNA', ...
        'N10 SNA', ...
        'N2 SNA'};

    na_patterns = { ...
        {{'NB','SNA'}, {}}, ...
        {{'N10','SNA'}, {}}, ...
        {{'N2','SNA'}, {}}};

    % ---------- FIGURE (2×1) ----------
    fig = figure('Color','w', ...
        'Name','RG over RB ratio (R8L+NEW, SNA+ only) – 2×1');
    tl = tiledlayout(fig, 2, 1, 'TileSpacing','compact', 'Padding','compact');

    % EndoH
    ax1 = nexttile(tl, 1);
    plotRatioRow(ax1, endo_labels, endo_patterns, keys, RG, RB, ...
                 ratioTrue, ratioPlot, ymax);
    title(ax1, 'EndoH series (SNA+)');
    box(ax1,'on');

    % NA
    ax2 = nexttile(tl, 2);
    plotRatioRow(ax2, na_labels, na_patterns, keys, RG, RB, ...
                 ratioTrue, ratioPlot, ymax);
    title(ax2, 'NA series (SNA+)');
    box(ax2,'on');
end

% =====================================================================
function plotRatioRow(ax, row_labels, row_patterns, keys, RG, RB, ...
                      ratioTrue, ratioPlot, ymax)

    n = numel(row_labels);
    x = 1:n;

    axes(ax); %#ok<LAXES>
    hold(ax,'on');

    for j = 1:n
        incl = row_patterns{j}{1};
        excl = row_patterns{j}{2};
        idx  = findKeyIdxByPattern(keys, incl, excl);

        if isnan(idx)
            y   = 0;
            txt = 'RG 0 / RB 0';
        else
            y   = ratioPlot(idx);
            RGc = RG(idx);
            RBc = RB(idx);

            if RBc > 0
                r = ratioTrue(idx);
                if isnan(r), r = 0; end
                txt = sprintf('RG %d / RB %d = %.2f', RGc, RBc, r);
            else
                if RGc > 0
                    txt = sprintf('RG %d / RB 0', RGc);
                else
                    txt = 'RG 0 / RB 0';
                end
            end
        end

        % Draw bar
        bar(ax, x(j), y, 0.6, 'FaceColor', [0.7 0.7 0.95], ...
                               'EdgeColor', [0.5 0.5 0.7]);

        % Text position
        if y == 0
            yText = 0.05 * ymax;
        else
            yText = y;
        end
        text(ax, x(j), yText, txt, ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom', ...
            'FontSize',7);
    end

    hold(ax,'off');

    ylim(ax, [0 ymax]);
    xlim(ax, [0.5 n+0.5]);
    xticks(ax, 1:n);
    xticklabels(ax, row_labels);
    xtickangle(ax, 45);
    ylabel(ax, 'RG (incl RGB) / RB-only');
end

% =====================================================================
function idx = findKeyIdxByPattern(keys, inclSubstrings, exclSubstrings)
    if nargin < 3 || isempty(exclSubstrings)
        exclSubstrings = {};
    end

    mask = true(size(keys));

    for i = 1:numel(inclSubstrings)
        s = inclSubstrings{i};
        mask = mask & contains(keys, s, 'IgnoreCase', true);
    end

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
