function plotColocCategoryBars_R8L_SNAplus(eventsRoot)
%PLOTCOLORCATEGORYBARS_R8L_SNAPLUS
%   3×1 stacked bar figure of coloc category fractions:
%       NoColoc / ColocGreen / ColocBlue / ColocAll
%
%   Only SNA+ samples, excluding R-L masters.
%   Conditions are defined via conditionKey(), pooling:
%       • E2 SNA  -> R8L 37C EndoH+ SNA+
%       • EB SNA  -> R8L 37C EndoH- SNA+

    if nargin < 1 || isempty(eventsRoot)
        eventsRoot = 'EventsAnalysis';
    end

    if ~isfolder(eventsRoot)
        error('Events root "%s" not found.', eventsRoot);
    end

    % -------- collect event files --------
    [eventFiles, condKeys, meta] = collectEventFiles(eventsRoot);

    % SNA+ only, non–R-L
    isSNA = contains({meta.sample}', 'SNA', 'IgnoreCase', true);
    isRL  = strcmp({meta.master}', 'R-L');
    keep  = isSNA & ~isRL;

    eventFiles = eventFiles(keep);
    condKeys   = condKeys(keep);
    meta       = meta(keep); %#ok<NASGU>

    if isempty(eventFiles)
        warning('No SNA+ non–R-L event files found for coloc category plot.');
        return;
    end

    % unique condition keys (after pooling)
    keys = unique(condKeys, 'stable');
    nK   = numel(keys);

    % fractions per key: [NoColoc, ColocGreen, ColocBlue, ColocAll]
    fractions = zeros(nK, 4);

    for i = 1:numel(eventFiles)
        k = find(strcmp(keys, condKeys{i}), 1);
        T = readtable(eventFiles{i});
        if isempty(T), continue; end

        CG = logical(T.ColocGreen);
        CB = logical(T.ColocBlue);
        CA = logical(T.ColocAll);
        NC = logical(T.NoColoc);

        fractions(k,1) = fractions(k,1) + nnz(NC);
        fractions(k,2) = fractions(k,2) + nnz(CG);
        fractions(k,3) = fractions(k,3) + nnz(CB);
        fractions(k,4) = fractions(k,4) + nnz(CA);
    end

    % convert to percentages
    rowSums = sum(fractions, 2);
    rowSums(rowSums == 0) = 1;  % avoid /0
    fractionsPct = 100 * fractions ./ rowSums;

    % ---------- row definitions (same as RG/RB plot) ----------

    % Row 1: temperature series
    row1_labels = { ...
        'R8L 25C EndoH- SNA+', ...
        'R8L 25C EndoH+ SNA+', ...
        'R8L 37C EndoH- SNA+', ...
        'R8L 37C EndoH+ SNA+', ...
        'R8L SNA+'};

    row1_patterns = { ...
        {{'R8L','25','EndoH-','SNA+'}, {}}, ...
        {{'R8L','25','EndoH+','SNA+'}, {}}, ...
        {{'R8L','37','EndoH-','SNA+'}, {}}, ...
        {{'R8L','37','EndoH+','SNA+'}, {}}, ...
        {{'R8L','SNA+'}, {'25','37','EndoH','NB','N10','N2','E10','EB','E2'}}};

    % Row 2: EndoH / glycosidase concentration
    row2_labels = { ...
        'R8L 37C EndoH- SNA+', ...
        'E10 SNA', ...
        'R8L 37C EndoH+ SNA+', ...
        'R8L SNA+'};

    row2_patterns = { ...
        {{'R8L','37','EndoH-','SNA+'}, {}}, ...
        {{'E10','SNA'}, {}}, ...
        {{'R8L','37','EndoH+','SNA+'}, {}}, ...
        {{'R8L','SNA+'}, {'25','37','EndoH','NB','N10','N2','E10','EB','E2'}}};

    % Row 3: NA concentration series
    row3_labels = { ...
        'NB SNA', ...
        'N10 SNA', ...
        'N2 SNA', ...
        'R8L SNA+'};

    row3_patterns = { ...
        {{'NB','SNA'}, {}}, ...
        {{'N10','SNA'}, {}}, ...
        {{'N2','SNA'}, {}}, ...
        {{'R8L','SNA+'}, {'25','37','EndoH','NB','N10','N2','E10','EB','E2'}}};

    % ---------- build figure ----------
    fig = figure('Color','w', ...
        'Name','Coloc categories per condition (R8L+NEW, SNA+ only) – 3×1');
    tl = tiledlayout(fig, 3, 1, 'TileSpacing','compact', 'Padding','compact');

    % Row 1
    ax1 = nexttile(tl, 1);
    rowFrac1 = buildRowFractions(row1_patterns, keys, fractionsPct);
    b1 = bar(ax1, rowFrac1, 'stacked');
    styleBars(b1);
    ylim(ax1, [0 100]);
    xlim(ax1, [0.5 size(rowFrac1,1)+0.5]);
    xticks(ax1, 1:size(rowFrac1,1));
    xticklabels(ax1, row1_labels);
    xtickangle(ax1, 45);
    ylabel(ax1, '% of events');
    title(ax1, 'Row 1 – Temperature series (R8L, SNA+)');
    box(ax1,'on');

    % Legend only on row 1
    legend(ax1, b1, {'NoColoc','ColocGreen','ColocBlue','ColocAll'}, ...
        'Location','northoutside','Orientation','horizontal');

    % Row 2
    ax2 = nexttile(tl, 2);
    rowFrac2 = buildRowFractions(row2_patterns, keys, fractionsPct);
    b2 = bar(ax2, rowFrac2, 'stacked');
    styleBars(b2);
    ylim(ax2, [0 100]);
    xlim(ax2, [0.5 size(rowFrac2,1)+0.5]);
    xticks(ax2, 1:size(rowFrac2,1));
    xticklabels(ax2, row2_labels);
    xtickangle(ax2, 45);
    ylabel(ax2, '% of events');
    title(ax2, 'Row 2 – EndoH / glycosidase concentration (R8L + E10, SNA+)');
    box(ax2,'on');

    % Row 3
    ax3 = nexttile(tl, 3);
    rowFrac3 = buildRowFractions(row3_patterns, keys, fractionsPct);
    b3 = bar(ax3, rowFrac3, 'stacked');
    styleBars(b3);
    ylim(ax3, [0 100]);
    xlim(ax3, [0.5 size(rowFrac3,1)+0.5]);
    xticks(ax3, 1:size(rowFrac3,1));
    xticklabels(ax3, row3_labels);
    xtickangle(ax3, 45);
    ylabel(ax3, '% of events');
    title(ax3, 'Row 3 – NA concentration series (NB/N10/N2, SNA+)');
    box(ax3,'on');

end

% =====================================================================
function rowFrac = buildRowFractions(row_patterns, keys, fractionsPct)
% Build an n×4 matrix of fractions for one row, in the desired order.

    n = numel(row_patterns);
    rowFrac = zeros(n, 4);

    for j = 1:n
        incl = row_patterns{j}{1};
        excl = row_patterns{j}{2};
        idx  = findKeyIdxByPattern(keys, incl, excl);

        if ~isnan(idx)
            rowFrac(j,:) = fractionsPct(idx,:);
        else
            rowFrac(j,:) = 0;
        end
    end
end

% =====================================================================
function idx = findKeyIdxByPattern(keys, inclSubstrings, exclSubstrings)
%FINDKEYIDXBYPATTERN  Find index in KEYS that matches all INCL and none of EXCL.

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

% =====================================================================
function styleBars(b)
%STYLEBARS  Apply consistent colors to stacked bar objects.

    if numel(b) < 4
        return;
    end

    b(1).FaceColor = [1 0 0];       % NoColoc – red
    b(2).FaceColor = [0 0.7 0];     % ColocGreen – green
    b(3).FaceColor = [0 0.4 1];     % ColocBlue – blue
    b(4).FaceColor = [0.6 0 0.6];   % ColocAll – purple
end
