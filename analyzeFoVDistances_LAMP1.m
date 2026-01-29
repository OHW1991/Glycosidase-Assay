function analyzeFoVDistances_LAMP1(rawDataRoot)
% Walk RawData and compute CoG distances for each FoV.
%
% See previous description for details.

    if nargin < 1 || isempty(rawDataRoot)
        rawDataRoot = 'RawData';
    end
    if ~isfolder(rawDataRoot)
        error('RawData root "%s" does not exist.', rawDataRoot);
    end

    % Thresholds
    redThresh_noSNA   = 3000;
    redThresh_withSNA = 4500;
    blueThresh        = 2000;
    minEventSize      = 8;   % pixels

    % Fresh AnalysisDistance
    analysisRoot = fullfile(pwd, 'AnalysisDistance');
    if isfolder(analysisRoot)
        try
            rmdir(analysisRoot,'s');
        catch ME
            warning('Could not remove existing AnalysisDistance: %s', ME.message);
        end
    end
    mkdir(analysisRoot);

    % ----- loop over date folders -----
    dateDirs = dir(rawDataRoot);
    dateDirs = dateDirs([dateDirs.isdir]);
    dateDirs = dateDirs(~ismember({dateDirs.name},{'.','..'}));

    for iD = 1:numel(dateDirs)
        dateName = dateDirs(iD).name;
        datePath = fullfile(rawDataRoot, dateName);

        lvl1 = dir(datePath);
        lvl1 = lvl1([lvl1.isdir]);
        lvl1 = lvl1(~ismember({lvl1.name},{'.','..'}));

        for i1 = 1:numel(lvl1)
            lvl1Name = lvl1(i1).name;
            lvl1Path = fullfile(datePath, lvl1Name);

            % >>> FIX: test lvl1Path itself, not its child <<<
            if isSampleFolder(lvl1Path)
                % lvl1 is already a sample folder
                processSampleFolder(analysisRoot, lvl1Path, lvl1Name, ...
                    redThresh_noSNA, redThresh_withSNA, blueThresh, minEventSize);
            else
                % lvl1 is a master (e.g. R8L, R-L); its children are samples
                sampleDirs = dir(lvl1Path);
                sampleDirs = sampleDirs([sampleDirs.isdir]);
                sampleDirs = sampleDirs(~ismember({sampleDirs.name},{'.','..'}));

                for iS = 1:numel(sampleDirs)
                    sampNameRaw = sampleDirs(iS).name;
                    sampPath    = fullfile(lvl1Path, sampNameRaw);
                    combinedName = sprintf('%s %s', lvl1Name, sampNameRaw);
                    processSampleFolder(analysisRoot, sampPath, combinedName, ...
                        redThresh_noSNA, redThresh_withSNA, blueThresh, minEventSize);
                end
            end
        end
    end

    fprintf('Done. Distance Excel files written under "%s".\n', analysisRoot);
end

% =======================================================================
% Does this folder itself look like a sample folder (contains FoVs)?
% =======================================================================
function tf = isSampleFolder(folderPath)
    tf = false;
    if ~isfolder(folderPath), return; end

    d = dir(folderPath);
    d = d([d.isdir]);
    d = d(~ismember({d.name},{'.','..'}));

    if isempty(d)
        % maybe TIFFs directly in this folder (we'll still treat it as sample)
        tifs = dir(fullfile(folderPath,'*.tif'));
        tf = ~isempty(tifs);
        return;
    end

    % Sample folders should have subdirs named like "FoV*"
    for i = 1:numel(d)
        if contains(d(i).name,'FoV','IgnoreCase',true)
            tf = true;
            return;
        end
    end

    % otherwise: assume not a sample
end

% =======================================================================
% Process a single sample folder (with FoV subfolders)
% =======================================================================
function processSampleFolder(analysisRoot, samplePath, sampleName, ...
    redThresh_noSNA, redThresh_withSNA, blueThresh, minEventSize)

    if ~isfolder(samplePath), return; end

    % Decide SNA+/SNA- from sampleName
    upperName   = upper(sampleName);
    hasSNAword  = contains(upperName,'SNA');
    hasSNAminus = contains(upperName,'SNA-');
    hasSNAplus  = contains(upperName,'SNA+');

    if hasSNAplus || (hasSNAword && ~hasSNAminus)
        redThresh = redThresh_withSNA;
    else
        redThresh = redThresh_noSNA;
    end

    outSampleDir = fullfile(analysisRoot, sampleName);
    if ~isfolder(outSampleDir)
        mkdir(outSampleDir);
    end

    % FoV dirs
    fovDirs = dir(samplePath);
    fovDirs = fovDirs([fovDirs.isdir]);
    fovDirs = fovDirs(~ismember({fovDirs.name},{'.','..'}));

    for iF = 1:numel(fovDirs)
        fovName = fovDirs(iF).name;
        fovPath = fullfile(samplePath, fovName);

        % skip non-FoV folders (just in case)
        if ~contains(fovName,'FoV','IgnoreCase',true)
            continue;
        end

        redFile  = findChannelTiff(fovPath, 1);
        blueFile = findChannelTiff(fovPath, 3);

        if isempty(redFile)
            fprintf('No red TIFF found in %s, skipping.\n', fovPath);
            continue;
        end

        R = imread(fullfile(fovPath, redFile));
        R = squeezeFirstPlane(R);
        if isempty(blueFile)
            B = zeros(size(R),'like',R);
        else
            B = imread(fullfile(fovPath, blueFile));
            B = squeezeFirstPlane(B);
        end

        R = double(R);
        B = double(B);

        redMask  = R >= redThresh;
        blueMask = B >= blueThresh;

        cc     = bwconncomp(redMask, 8);
        imSize = size(R);

        % *** use consistent names here ***
        eventID = [];
        redX    = [];
        redY    = [];
        rbX     = [];
        rbY     = [];
        distPix = [];

        evCount = 0;
        for ev = 1:cc.NumObjects
            idx = cc.PixelIdxList{ev};
            if numel(idx) < minEventSize
                continue;
            end

            % -----------------------------
            % Intensity-weighted RED CoG
            % -----------------------------
            [rows, cols] = ind2sub(imSize, idx);
            Ired = R(idx);              % already double
            wR   = max(Ired, eps);      % avoid all-zero weights
            sumW = sum(wR);

            CoG_red_x = sum(cols .* wR) / sumW;
            CoG_red_y = sum(rows .* wR) / sumW;

            % -----------------------------
            % Intensity-weighted BLUE CoG on RB pixels
            % -----------------------------
            rbMaskLocal = blueMask(idx);
            rbIdx       = idx(rbMaskLocal);

            if isempty(rbIdx)
                CoG_rb_x = NaN;
                CoG_rb_y = NaN;
                distRB   = NaN;
            else
                [rbRows, rbCols] = ind2sub(imSize, rbIdx);
                Iblue = B(rbIdx);
                wB    = max(Iblue, eps);
                sumWB = sum(wB);

                CoG_rb_x = sum(rbCols .* wB) / sumWB;
                CoG_rb_y = sum(rbRows .* wB) / sumWB;

                dx = CoG_rb_x - CoG_red_x;
                dy = CoG_rb_y - CoG_red_y;
                distRB = sqrt(dx.^2 + dy.^2);
            end

            % Store
            evCount              = evCount + 1;
            eventID(evCount)     = ev;         %#ok<AGROW>
            redX(evCount)        = CoG_red_x;  %#ok<AGROW>
            redY(evCount)        = CoG_red_y;  %#ok<AGROW>
            rbX(evCount)         = CoG_rb_x;   %#ok<AGROW>
            rbY(evCount)         = CoG_rb_y;   %#ok<AGROW>
            distPix(evCount)     = distRB;     %#ok<AGROW>
        end

        % ---- build table ----
        if isempty(eventID)
            T = table([],[],[],[],[],[], ...
                'VariableNames',{'EventID','RedX','RedY','RB_X','RB_Y','DistancePixels'});
        else
            T = table(eventID(:), redX(:), redY(:), rbX(:), rbY(:), distPix(:), ...
                'VariableNames',{'EventID','RedX','RedY','RB_X','RB_Y','DistancePixels'});
        end

        outFile = fullfile(outSampleDir, [fovName '_distances.xlsx']);
        try
            if isfile(outFile), delete(outFile); end
            writetable(T,outFile);
        catch ME
            warning('Failed writing %s: %s', outFile, ME.message);
        end
    end
end

% =======================================================================
% Find TIFF file for a given channel index in a FoV folder
% =======================================================================
function fname = findChannelTiff(fovPath, chanIdx)
    fname = '';

    % Common patterns: *_c1*.tif, *_Channel1*.tif
    pat1 = sprintf('*c%d*.tif',chanIdx);
    files = dir(fullfile(fovPath,pat1));
    if ~isempty(files)
        fname = files(1).name;
        return;
    end

    pat2 = sprintf('*Channel%d*.tif',chanIdx);
    files = dir(fullfile(fovPath,pat2));
    if ~isempty(files)
        fname = files(1).name;
        return;
    end

    % Fallback: just take sorted TIFFs and map:
    %   first = red, last = blue
    allTifs = dir(fullfile(fovPath,'*.tif'));
    if isempty(allTifs)
        return;
    end
    [~,idx] = sort({allTifs.name});
    allTifs = allTifs(idx);

    if chanIdx == 1
        fname = allTifs(1).name;
    elseif chanIdx == 3
        fname = allTifs(end).name;
    else
        fname = '';
    end
end

% =======================================================================
% Max-projection if stack
% =======================================================================
function Iout = squeezeFirstPlane(I)
    if ndims(I) > 2
        Iout = max(I,[],3);
    else
        Iout = I;
    end
end
