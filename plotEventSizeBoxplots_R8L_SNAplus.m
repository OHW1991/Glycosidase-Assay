function plotEventSizeBoxplots_R8L_SNAplus(eventsRoot)
%PLOTEVENTSIZEBOXPLOTS_R8L_SNAPLUS
%   Event size (NumPixels) per conditionKey,
%   only SNA+ and non–R-L.

if nargin < 1 || isempty(eventsRoot)
    eventsRoot = 'EventsAnalysis';
end

[eventFiles, condKeys, meta] = collectEventFiles(eventsRoot);

isSNA = contains({meta.sample}','SNA','IgnoreCase',true);
isRL  = strcmp({meta.master}','R-L');
keep  = isSNA & ~isRL;

eventFiles = eventFiles(keep);
condKeys   = condKeys(keep);

keys = unique(condKeys,'stable');
nK   = numel(keys);
sizes = cell(nK,1);

for i = 1:numel(eventFiles)
    k = find(strcmp(keys,condKeys{i}),1);
    T = readtable(eventFiles{i});
    if isempty(T), continue; end
    sizes{k} = [sizes{k}; double(T.NumPixels)];
end

figure('Color','w','Name','Event size (R8L+NEW, SNA+ only)');
hold on;
for k = 1:nK
    if isempty(sizes{k}), continue; end
    boxchart(k*ones(size(sizes{k})), sizes{k}, 'BoxWidth',0.5);
end
hold off;
box on;

xlim([0.5 nK+0.5]);
xticks(1:nK); xticklabels(keys); xtickangle(45);
ylabel('Event size (pixels)');
title('Red event size distribution (R8L+NEW, SNA+ only)');

end
