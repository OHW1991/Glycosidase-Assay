function plotEventSizeBoxplots(eventsRoot)
%PLOTEVENTSIZEBOXPLOTS  Boxplots of event sizes (NumPixels) per condition.

if nargin < 1 || isempty(eventsRoot)
    eventsRoot = 'EventsAnalysis';
end

[eventFiles, condKeys] = collectEventFiles(eventsRoot);

keys = unique(condKeys,'stable');
nK   = numel(keys);
sizes = cell(nK,1);

for i = 1:numel(eventFiles)
    k = find(strcmp(keys,condKeys{i}),1);
    T = readtable(eventFiles{i});
    if isempty(T), continue; end
    sizes{k} = [sizes{k}; double(T.NumPixels)];
end

figure('Color','w','Name','Event size per condition');
hold on;
for k = 1:nK
    if isempty(sizes{k}), continue; end
    boxchart(k*ones(size(sizes{k})), sizes{k}, 'BoxWidth',0.5);
end
hold off;

xlim([0.5 nK+0.5]);
xticks(1:nK); xticklabels(keys); xtickangle(45);
ylabel('Event size (pixels)');
title('Red event size distribution');

end
