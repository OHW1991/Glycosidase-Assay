function plotRGandRB_counts(eventsRoot)
%PLOTRGANDRB_COUNTS  Plot RG & RB-only counts per conditionKey.
%
%   plotRGandRB_counts('EventsAnalysis');

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

allCounts = [RG,RB];
maxCount  = max(allCounts);
if isempty(maxCount) || maxCount<=0, maxCount = 1; end

figure('Color','w','Name','RG and RB counts per condition');
hold on;

edgeCol   = [0.4 0 0];
faceGreen = [0 0.7 0];
faceBlue  = [0 0.4 1];

x = 1:nK;
for k = 1:nK
    % RG
    plot(x(k), RG(k), 'o', 'MarkerFaceColor',faceGreen, ...
        'MarkerEdgeColor',edgeCol, 'MarkerSize',8,'LineWidth',1.5);
    % RB
    plot(x(k), RB(k), 'o', 'MarkerFaceColor',faceBlue, ...
        'MarkerEdgeColor',edgeCol, 'MarkerSize',8,'LineWidth',1.5);
end

ylim([0 maxCount*1.2]);
xlim([0.5 nK+0.5]);
xticks(1:nK); xticklabels(keys); xtickangle(45);
ylabel('Number of events');
title('RG (incl RGB) and RB-only counts');

hRG = plot(nan,nan,'o','MarkerFaceColor',faceGreen, ...
    'MarkerEdgeColor',edgeCol,'MarkerSize',8,'LineWidth',1.5);
hRB = plot(nan,nan,'o','MarkerFaceColor',faceBlue, ...
    'MarkerEdgeColor',edgeCol,'MarkerSize',8,'LineWidth',1.5);
legend([hRG,hRB],{'RG (incl RGB)','RB-only'}, ...
    'Location','northoutside','Orientation','horizontal');

hold off;

end
