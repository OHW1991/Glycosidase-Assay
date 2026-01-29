function analyzeFoVEvents(map, eventsRoot)
%ANALYZEFOVEVENTS  Detect red events per FoV and write *_events.xlsx tables.
%
%   analyzeFoVEvents(map, 'EventsAnalysis');
%
% Uses different thresholds for 20251125 (old) vs 20251130 (new):
%
%   OLD (20251125):
%       redThr (SNA-) = 2000
%       redThr (SNA+) = 4500
%       greenThr      = 8000
%       blueThr       = 2000
%
%   NEW (20251130):
%       redThr (SNA-) = 3000
%       redThr (SNA+) = 6000
%       greenThr      = 4000
%       blueThr       = 2000
%
% Event = connected component of >= 8 red pixels above redThr (8-connectivity).
%
% For each event we compute:
%   NumPixels, NumGreen, NumBlue, NumBoth,
%   ColocGreen, ColocBlue, ColocAll, NoColoc
%
% RG (red-green)   = pixels where red>thr & green>thr
% RB (red-blue)    = pixels where red>thr & blue>thr
% Both             = pixels with red>thr & green>thr & blue>thr
%
% ColocGreen = (NumGreen / NumPixels) >= 0.75
% ColocBlue  = (NumBlue  / NumPixels) >= 0.15
% ColocAll   = ColocGreen & ColocBlue
% NoColoc    = ~ColocGreen & ~ColocBlue

    if nargin < 2 || isempty(eventsRoot)
        eventsRoot = 'EventsAnalysis';
    end
    if ~isfolder(eventsRoot)
        mkdir(eventsRoot);
    end

    for i = 1:numel(map)
        rec = map(i);

        % ----- thresholds by date / SNA -----
        switch rec.date
            case '20251125'
                if rec.isSNA
                    redThr = 4500;
                else
                    redThr = 3000; %2000
                end
                greenThr = 8000;
            case '20251130'
                if rec.isSNA
                    redThr = 4500; %6000
                else
                    redThr = 3000;
                end
                greenThr = 4000;
            otherwise
                warning('Unknown date %s, skipping', rec.date);
                continue;
        end
        blueThr = 2000;
        sizeThr = 8;

        % ----- output folder -----
        if isempty(rec.master)
            condPath = fullfile(eventsRoot, rec.date, rec.sample);
        else
            condPath = fullfile(eventsRoot, rec.date, rec.master, rec.sample);
        end
        if ~isfolder(condPath)
            mkdir(condPath);
        end

        % ----- read channels -----
        R = double(imread(rec.c1));
        G = double(imread(rec.c2));
        B = double(imread(rec.c3));

        redMask = R > redThr;
        CC = bwconncomp(redMask, 8);

        rows = [];   % will hold event summary rows

        greenMask = G > greenThr;
        blueMask  = B > blueThr;

        for e = 1:CC.NumObjects
            pixIdx = CC.PixelIdxList{e};
            if numel(pixIdx) < sizeThr
                continue;
            end

            eventMask = false(size(redMask));
            eventMask(pixIdx) = true;

            numPixels = numel(pixIdx);
            numGreen  = nnz(eventMask & greenMask);
            numBlue   = nnz(eventMask & blueMask);
            numBoth   = nnz(eventMask & greenMask & blueMask);

            pG = numGreen / numPixels;
            pB = numBlue  / numPixels;

            colocGreen = pG >= 0.75;
            colocBlue  = pB >= 0.15;
            colocAll   = colocGreen && colocBlue;
            noColoc    = ~colocGreen && ~colocBlue;

            rows = [rows; ...
                e, numPixels, numGreen, numBlue, numBoth, ...
                colocGreen, colocBlue, colocAll, noColoc]; %#ok<AGROW>
        end

        % ----- handle FoVs with ZERO events -----
        if isempty(rows)
            % make it a 0-by-9 matrix so array2table knows there are 9 vars
            rows = zeros(0,9);
        end

        T = array2table(rows, ...
            'VariableNames', {'EventID','NumPixels','NumGreen','NumBlue','NumBoth', ...
                              'ColocGreen','ColocBlue','ColocAll','NoColoc'});

        outName = sprintf('%s_%s_events.xlsx', rec.sample, rec.fov);
        outPath = fullfile(condPath, outName);
        writetable(T, outPath);
    end
end
