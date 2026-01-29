function KW_greenProportion_R8L_SNAplus(eventsRoot)
%KW_GREENPROPORTION_R8L_SNAPLUS
%   Kruskal–Wallis on per-FoV green proportion (RG / (RG+RB)),
%   only SNA+ and not R-L.

if nargin < 1 || isempty(eventsRoot)
    eventsRoot = 'EventsAnalysis';
end

[eventFiles, condKeys, meta] = collectEventFiles(eventsRoot);

isSNA = contains({meta.sample}','SNA','IgnoreCase',true);
isRL  = strcmp({meta.master}','R-L');
keep  = isSNA & ~isRL;

eventFiles = eventFiles(keep);
condKeys   = condKeys(keep);
meta       = meta(keep); %#ok<NASGU>  % kept for symmetry, not used

keys = unique(condKeys,'stable');
nK   = numel(keys);

all_pG    = [];
all_group = [];

for i = 1:numel(eventFiles)
    k = find(strcmp(keys,condKeys{i}),1);
    T = readtable(eventFiles{i});
    if isempty(T), continue; end

    CG = logical(T.ColocGreen);
    CB = logical(T.ColocBlue);
    CA = logical(T.ColocAll);

    RG = nnz(CG | CA);
    RB = nnz(CB & ~CA);

    if RG + RB == 0
        continue;
    end

    pG = RG / (RG + RB);
    all_pG(end+1)    = pG; %#ok<AGROW>
    all_group(end+1) = k; %#ok<AGROW>
end

if isempty(all_pG)
    warning('No usable FoV datapoints for KW (R8L+NEW, SNA+ only).');
    return;
end

[pKW, ~, stats] = kruskalwallis(all_pG, all_group, 'off');
fprintf('KW p-value (R8L+NEW, SNA+ only) = %.4g\n', pKW);

c = multcompare(stats, 'CType','dunn-sidak','Display','off');
disp(array2table(c, ...
    'VariableNames',{'Group1','Group2','LowerCI','Diff','UpperCI','pValue'}));

figure('Color','w','Name','KW green proportion (R8L+NEW, SNA+ only)');
boxchart(all_group, all_pG, 'BoxFaceAlpha',0.4);
hold on;
swarmchart(all_group, all_pG, 'filled');
hold off;
box on;

xticks(1:nK); xticklabels(keys); xtickangle(45);
ylabel('Proportion green (RG / (RG + RB))');
ylim([0 1]);
title('Kruskal–Wallis on green proportion (R8L+NEW, SNA+ only)');

end
