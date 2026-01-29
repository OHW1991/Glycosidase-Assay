function plotColocFractionBoxplots_R8L_SNAplus(eventsRoot)
%PLOTCOLOCFRACTIONBOXPLOTS_R8L_SNAPLUS
%   Like plotColocFractionBoxplots, but:
%       • only SNA+ samples
%       • excludes any R-L master
%       • pools E2 SNA with R8L 37C EndoH+ SNA+
%         and EB SNA with R8L 37C EndoH- SNA+ (via conditionKey)
%
%   Fractions are per-event % of red pixels in:
%       Red-only, RB-only, RG-only, Triple

if nargin < 1 || isempty(eventsRoot)
    eventsRoot = 'EventsAnalysis';
end

[eventFiles, condKeys, meta] = collectEventFiles(eventsRoot);

% filter to SNA+ and not R-L
isSNA = contains({meta.sample}','SNA','IgnoreCase',true);
isRL  = strcmp({meta.master}','R-L');
keep  = isSNA & ~isRL;

eventFiles = eventFiles(keep);
condKeys   = condKeys(keep);

keys = unique(condKeys,'stable');
nK   = numel(keys);

redOnlyAll  = cell(nK,1);
redBlueAll  = cell(nK,1);
redGreenAll = cell(nK,1);
tripleAll   = cell(nK,1);

for i = 1:numel(eventFiles)
    key = condKeys{i};
    k   = find(strcmp(keys,key),1);

    T = readtable(eventFiles{i});
    if isempty(T), continue; end

    N    = double(T.NumPixels);
    G    = double(T.NumGreen);
    B    = double(T.NumBlue);
    Both = double(T.NumBoth);

    RG_only = max(G - Both,0);
    RB_only = max(B - Both,0);
    Triple  = Both;
    RedOnly = max(N - RG_only - RB_only - Triple,0);

    valid = N > 0;
    N      = N(valid);
    RedOnly= RedOnly(valid);
    RG_only= RG_only(valid);
    RB_only= RB_only(valid);
    Triple = Triple(valid);

    redOnlyAll{k}  = [redOnlyAll{k};  100*(RedOnly./N)];
    redBlueAll{k}  = [redBlueAll{k};  100*(RB_only./N)];
    redGreenAll{k} = [redGreenAll{k}; 100*(RG_only./N)];
    tripleAll{k}   = [tripleAll{k};   100*(Triple./N)];
end

figure('Color','w','Name','Coloc fractions (R8L+NEW, SNA+ only)');
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

% 1) Red-only
nexttile;
hold on;
for k = 1:nK
    if isempty(redOnlyAll{k}), continue; end
    boxchart(k*ones(size(redOnlyAll{k})), redOnlyAll{k}, 'BoxWidth',0.5);
end
hold off;
box on;
xlim([0.5 nK+0.5]); ylim([0 100]);
xticks(1:nK); xticklabels(keys); xtickangle(45);
ylabel('% red-only'); title('Red-only');

% 2) RB-only
nexttile;
hold on;
for k = 1:nK
    if isempty(redBlueAll{k}), continue; end
    boxchart(k*ones(size(redBlueAll{k})), redBlueAll{k}, 'BoxWidth',0.5);
end
hold off;
box on;
xlim([0.5 nK+0.5]); ylim([0 100]);
xticks(1:nK); xticklabels(keys); xtickangle(45);
ylabel('% RB-only'); title('Red-Blue only');

% 3) RG-only
nexttile;
hold on;
for k = 1:nK
    if isempty(redGreenAll{k}), continue; end
    boxchart(k*ones(size(redGreenAll{k})), redGreenAll{k}, 'BoxWidth',0.5);
end
hold off;
box on;
xlim([0.5 nK+0.5]); ylim([0 100]);
xticks(1:nK); xticklabels(keys); xtickangle(45);
ylabel('% RG-only'); title('Red-Green only');

% 4) Triple
nexttile;
hold on;
for k = 1:nK
    if isempty(tripleAll{k}), continue; end
    boxchart(k*ones(size(tripleAll{k})), tripleAll{k}, 'BoxWidth',0.5);
end
hold off;
box on;
xlim([0.5 nK+0.5]); ylim([0 100]);
xticks(1:nK); xticklabels(keys); xtickangle(45);
ylabel('% triple'); title('Triple label');

end
