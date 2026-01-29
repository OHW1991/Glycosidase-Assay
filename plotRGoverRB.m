function plotRGoverRB(eventsRoot)
%PLOTRGOVERRB  Plot RG(incl triple) / RB-only ratio per conditionKey.
%
%   plotRGoverRB('EventsAnalysis');
%
% Uses pooled counts (all FoVs per conditionKey).

if nargin < 1 || isempty(eventsRoot)
    eventsRoot = 'EventsAnalysis';
end

[eventFiles, condKeys] = collectEventFiles(eventsRoot);
keys = unique(condKeys,'stable');
nK   = numel(keys);

RG = zeros(1,nK);
RB = zeros(1,nK);

for i = 1:numel(eventFiles)
    k = find(strcmp(keys,condKeys{i}),1);
    T = readtable(eventFiles{i});
    if isempty(T), continue; end

    CG = logical(T.ColocGreen);
    CB = logical(T.ColocBlue);
    CA = logical(T.ColocAll);

    RB(k) = RB(k) + nnz(CB & ~CA);   % RB-only
    RG(k) = RG(k) + nnz(CG | CA);    % RG + triple
end

ratio = nan(1,nK);
hasRB = RB > 0;
ratio(hasRB) = RG(hasRB) ./ RB(hasRB);

% If RB==0 & RG>0, set bar slightly above max finite
maxFinite = max(ratio(hasRB));
if isempty(maxFinite) || isnan(maxFinite), maxFinite = 1; end
idxInf = (RB==0 & RG>0);
ratio(idxInf) = maxFinite * 1.05;

ymax = max(ratio(~isnan(ratio)));
if isempty(ymax) || ymax <= 0, ymax = 1; end
ymax = ymax * 1.2;

figure('Color','w','Name','RG/RB ratio per condition');
b = bar(ratio);
b.FaceColor = [0.6 0.6 0.9];
xlim([0.5 nK+0.5]);
ylim([0 ymax]);
xticks(1:nK); xticklabels(keys); xtickangle(45);
ylabel('RG (incl RGB) / RB-only');
title('RG over RB ratio per conditionKey');

% annotate with counts
hold on;
for k = 1:nK
    if isnan(ratio(k)), continue; end
    txt = sprintf('RG %d / RB %d', RG(k), RB(k));
    y   = ratio(k);
    if y == 0, y = 0.05*ymax; end
    text(k, y, txt, 'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom','FontSize',7);
end
hold off;

end
