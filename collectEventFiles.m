function [eventFiles, condKeys, meta] = collectEventFiles(eventsRoot)
%COLLECTEVENTFILES  Return list of *_events.xlsx and condition keys.
%
%   [files, keys, meta] = collectEventFiles('EventsAnalysis');
%
% Condition key pools:
%   - 20251125 / R8L / 'R8L 37C EndoH+ SNA+'  +  20251130 / 'E2 SNA'
%   - 20251125 / R8L / 'R8L 37C EndoH- SNA+'  +  20251130 / 'EB SNA'

if nargin < 1 || isempty(eventsRoot)
    eventsRoot = 'EventsAnalysis';
end

eventFiles = {};
condKeys   = {};
meta       = struct('date',{},'master',{},'sample',{});

dateDirs = dir(eventsRoot);
dateDirs = dateDirs([dateDirs.isdir]);
dateDirs = dateDirs(~ismember({dateDirs.name},{'.','..'}));

for iD = 1:numel(dateDirs)
    dateName = dateDirs(iD).name;
    datePath = fullfile(eventsRoot, dateName);

    switch dateName
        case '20251125'
            masterDirs = dir(datePath);
            masterDirs = masterDirs([masterDirs.isdir]);
            masterDirs = masterDirs(~ismember({masterDirs.name},{'.','..'}));

            for iM = 1:numel(masterDirs)
                masterName = masterDirs(iM).name;
                masterPath = fullfile(datePath, masterName);

                sampleDirs = dir(masterPath);
                sampleDirs = sampleDirs([sampleDirs.isdir]);
                sampleDirs = sampleDirs(~ismember({sampleDirs.name},{'.','..'}));

                for iS = 1:numel(sampleDirs)
                    sampleName = sampleDirs(iS).name;
                    samplePath = fullfile(masterPath, sampleName);

                    xlsx = dir(fullfile(samplePath,'*_events.xlsx'));
                    for iF = 1:numel(xlsx)
                        eventFiles{end+1,1} = fullfile(samplePath, xlsx(iF).name); %#ok<AGROW>
                        condKeys{end+1,1}   = conditionKey(dateName, masterName, sampleName);
                        meta(end+1) = struct('date',dateName,'master',masterName, ...
                                             'sample',sampleName); %#ok<AGROW>
                    end
                end
            end

        case '20251130'
            sampleDirs = dir(datePath);
            sampleDirs = sampleDirs([sampleDirs.isdir]);
            sampleDirs = sampleDirs(~ismember({sampleDirs.name},{'.','..'}));

            for iS = 1:numel(sampleDirs)
                sampleName = sampleDirs(iS).name;
                % Stock was already not mapped in mapRawData -> no events
                samplePath = fullfile(datePath, sampleName);

                if ~isfolder(samplePath), continue; end

                xlsx = dir(fullfile(samplePath,'*_events.xlsx'));
                for iF = 1:numel(xlsx)
                    eventFiles{end+1,1} = fullfile(samplePath, xlsx(iF).name); %#ok<AGROW>
                    condKeys{end+1,1}   = conditionKey(dateName, '', sampleName);
                    meta(end+1) = struct('date',dateName,'master','', ...
                                         'sample',sampleName); %#ok<AGROW>
                end
            end
    end
end

end

function key = conditionKey(dateName, masterName, sampleName)
% Canonical condition name with pooling of 37C EndoH± SNA+ and E2/EB.

% Normalize spaces
sampleName = strtrim(sampleName);

if strcmp(dateName,'20251125')
    key = sprintf('%s | %s', masterName, sampleName);

elseif strcmp(dateName,'20251130')
    % pool E2 SNA with old R8L 37C EndoH+ SNA+
    if strcmpi(sampleName,'E2 SNA')
        key = 'R8L | R8L 37C EndoH+ SNA+';
    elseif strcmpi(sampleName,'EB SNA')
        key = 'R8L | R8L 37C EndoH- SNA+';
    else
        key = sprintf('NEW | %s', sampleName);
    end
else
    key = sprintf('%s | %s | %s', dateName, masterName, sampleName);
end

end
