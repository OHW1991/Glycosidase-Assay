function map = mapRawData(rawRoot)
%MAPRAWDATA  Map all FoVs and channels in RawData.
%
%   map = mapRawData('RawData');
%
% map is a struct array with fields:
%   date      - '20251125' or '20251130'
%   master    - 'R8L', 'R-L', or '' (for 20251130 samples)
%   sample    - sample/condition folder name
%   fov       - FoV folder name
%   c1, c2, c3 - full paths to tiff images
%   isSNA     - logical, SNA+ if name contains 'SNA' or 'SNA+'

if nargin < 1 || isempty(rawRoot)
    rawRoot = 'RawData';
end

if ~isfolder(rawRoot)
    error('RawData folder "%s" not found.', rawRoot);
end

map = struct('date',{},'master',{},'sample',{},'fov',{}, ...
             'c1',{},'c2',{},'c3',{},'isSNA',{});

dateDirs = dir(rawRoot);
dateDirs = dateDirs([dateDirs.isdir]);
dateDirs = dateDirs(~ismember({dateDirs.name},{'.','..'}));

for iD = 1:numel(dateDirs)
    dateName = dateDirs(iD).name;
    datePath = fullfile(rawRoot, dateName);

    switch dateName
        case '20251125'  % old analysis: R8L and R-L masters
            masterDirs = dir(datePath);
            masterDirs = masterDirs([masterDirs.isdir]);
            masterDirs = masterDirs(~ismember({masterDirs.name},{'.','..'}));

            for iM = 1:numel(masterDirs)
                masterName = masterDirs(iM).name;   % 'R8L' or 'R-L'
                masterPath = fullfile(datePath, masterName);

                sampleDirs = dir(masterPath);
                sampleDirs = sampleDirs([sampleDirs.isdir]);
                sampleDirs = sampleDirs(~ismember({sampleDirs.name},{'.','..'}));

                for iS = 1:numel(sampleDirs)
                    sampleName = sampleDirs(iS).name;
                    samplePath = fullfile(masterPath, sampleName);

                    fovDirs = dir(samplePath);
                    fovDirs = fovDirs([fovDirs.isdir]);
                    fovDirs = fovDirs(~ismember({fovDirs.name},{'.','..'}));

                    for iF = 1:numel(fovDirs)
                        fovName = fovDirs(iF).name;
                        fovPath = fullfile(samplePath, fovName);

                        c1 = findTiff(fovPath,'c1');
                        c2 = findTiff(fovPath,'c2');
                        c3 = findTiff(fovPath,'c3');

                        if isempty(c1) || isempty(c2) || isempty(c3)
                            warning('Missing channels in %s', fovPath);
                            continue;
                        end

                        isSNA = contains(sampleName,'SNA','IgnoreCase',true);

                        map(end+1) = struct( ...
                            'date',   dateName, ...
                            'master', masterName, ...
                            'sample', sampleName, ...
                            'fov',    fovName, ...
                            'c1',     c1, ...
                            'c2',     c2, ...
                            'c3',     c3, ...
                            'isSNA',  isSNA); %#ok<AGROW>
                    end
                end
            end

        case '20251130'  % new analysis: 14 sample folders
            sampleDirs = dir(datePath);
            sampleDirs = sampleDirs([sampleDirs.isdir]);
            sampleDirs = sampleDirs(~ismember({sampleDirs.name},{'.','..'}));

            for iS = 1:numel(sampleDirs)
                sampleName = sampleDirs(iS).name;

                % ignore Stock and Stock SNA completely
                if strcmpi(sampleName,'Stock') || strcmpi(sampleName,'Stock SNA')
                    continue;
                end

                samplePath = fullfile(datePath, sampleName);

                fovDirs = dir(samplePath);
                fovDirs = fovDirs([fovDirs.isdir]);
                fovDirs = fovDirs(~ismember({fovDirs.name},{'.','..'}));

                for iF = 1:numel(fovDirs)
                    fovName = fovDirs(iF).name;
                    fovPath = fullfile(samplePath, fovName);

                    c1 = findTiff(fovPath,'c1');
                    c2 = findTiff(fovPath,'c2');
                    c3 = findTiff(fovPath,'c3');

                    if isempty(c1) || isempty(c2) || isempty(c3)
                        warning('Missing channels in %s', fovPath);
                        continue;
                    end

                    isSNA = contains(sampleName,'SNA','IgnoreCase',true);

                    map(end+1) = struct( ...
                        'date',   dateName, ...
                        'master', '', ...
                        'sample', sampleName, ...
                        'fov',    fovName, ...
                        'c1',     c1, ...
                        'c2',     c2, ...
                        'c3',     c3, ...
                        'isSNA',  isSNA); %#ok<AGROW>
                end
            end
        otherwise
            warning('Unexpected date folder: %s (ignored)', dateName);
    end
end

end

function path = findTiff(folder, pattern)
    d = dir(fullfile(folder, sprintf('*%s*.tif',pattern)));
    if isempty(d)
        path = '';
    else
        path = fullfile(folder,d(1).name);
    end
end
