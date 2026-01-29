function plotColocCategoryBars(eventsRoot)
%PLOTCOLORCATEGORYBARS  Stacked bar of coloc categories per condition.
%
% Categories:
%   NoColoc   (red-only)
%   ColocGreen (>=75% green)
%   ColocBlue  (>=15% blue)
%   ColocAll   (both)
%
% The four bars are stacked to 100% for each condition.

if nargin < 1 || isempty(eventsRoot)
    eventsRoot = 'EventsAnalysis';
end

[eventFiles, condKeys] = collectEventFiles(eventsRoot);

keys = unique(condKeys,'stable');
nK   = numel(keys);

fractions = zeros(nK,4); % [NoColoc, ColocGreen, ColocBlue, ColocAll]

for i = 1:numel(eventFiles)
    k = find(strcmp(keys,condKeys{i}),1);
    T = readtable(eventFiles{i});
    if isempty(T), continue; end

    CG = logical(T.ColocGreen);
    CB = logical(T.ColocBlue);
    CA = logical(T.ColocAll);
    NC = logical(T.NoColoc);

    nEvents = height(T);
    if nEvents == 0, continue; end

    fractions(k,1) = fractions(k,1) + nnz(NC);
    fractions(k,2) = fractions(k,2) + nnz(CG);
    fractions(k,3) = fractions(k,3) + nnz(CB);
    fractions(k,4) = fractions(k,4) + nnz(CA);
end

rowSums = sum(fractions,2);
rowSums(rowSums==0) = 1;
fractionsPct = 100 * fractions ./ rowSums;

figure('Color','w','Name','Coloc category fractions');
b = bar(fractionsPct, 'stacked');
b(1).FaceColor = [1 0 0];       % NoColoc - red
b(2).FaceColor = [0 0.7 0];     % ColocGreen - green
b(3).FaceColor = [0 0.4 1];     % ColocBlue - blue
b(4).FaceColor = [0.6 0 0.6];   % ColocAll - purple

ylim([0 100]);
xlim([0.5 nK+0.5]);
xticks(1:nK); xticklabels(keys); xtickangle(45);
ylabel('% of events');
legend({'NoColoc','ColocGreen','ColocBlue','ColocAll'}, ...
    'Location','northoutside','Orientation','horizontal');
title('Coloc categories per condition');

end
