function KW_greenProportion(eventsRoot)
%KW_GREENPROPORTION  KW on green proportion per conditionKey (SNA+ only).
%
%   KW_greenProportion('EventsAnalysis');
%
% For each FoV:
%   RGplus = # (ColocGreen OR ColocAll)
%   RBonly = # (ColocBlue AND ~ColocAll)
%   pG     = RGplus / (RGplus + RBonly)
%
% Only SNA+ samples are used (sample or conditionKey contains 'SNA+').

if nargin < 1 || isempty(eventsRoot)
    eventsRoot = 'EventsAnalysis';
end

[eventFiles, condKeys, meta] = collectEventFiles(eventsRoot);

% restrict to SNA+ (by sample name, contained in meta.sample)
isSNAplus = contains({meta.sample}','SNA','IgnoreCase',true);
eventFiles = eventFiles(isSNAplus);
condKeys   = condKeys(isSNAplus);
meta       = meta(isSNAplus);

keys = unique(condKeys,'stable');
nK   = numel(keys);

all_pG   = [];
all_group = [];

for i = 1:numel(eventFiles)
    key = condKeys{i};
    k   = find(strcmp(keys,key),1);

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
    warning('No usable FoV datapoints for KW.');
    return;
end

[pKW, ~, stats] = kruskalwallis(all_pG, all_group, 'off');
fprintf('KW p-value (SNA+ only) = %.4g\n', pKW);

c = multcompare(stats, 'CType','dunn-sidak','Display','off');
disp(array2table(c, ...
    'VariableNames',{'Group1','Group2','LowerCI','Diff','UpperCI','pValue'}));

figure('Color','w','Name','KW green proportion (SNA+ only)');
boxchart(all_group, all_pG, 'BoxFaceAlpha',0.4);
hold on;
swarmchart(all_group, all_pG, 'filled');
hold off;

xticks(1:nK); xticklabels(keys); xtickangle(45);
ylabel('Proportion green (RG / (RG + RB))');
ylim([0 1]);
title('Kruskal–Wallis on green proportion (SNA+ only)');

end
