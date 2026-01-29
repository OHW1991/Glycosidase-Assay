function plotCoGDeltaTrend_LAMP1(analysisRoot)
%PLOTCOGDELTATREND_LAMP1
%   Third figure: trendlines of ΔD (mean SNA+ − mean SNA−) vs enzyme level
%   for:
%       • EndoH series at 37°C   (EndoH−, 1:20, 1:2)
%       • NA series              (NB, N10, N2)
%
%   X-ticks:  "No enzyme", "Low enzyme", "High enzyme".
%
%   Usage:
%       plotCoGDeltaTrend_LAMP1('AnalysisDistance');

    if nargin < 1 || isempty(analysisRoot)
        analysisRoot = 'AnalysisDistance';
    end
    if ~isfolder(analysisRoot)
        error('Analysis root "%s" not found.', analysisRoot);
    end

    % pixels → µm (same as other plots)
    pixToUm = 101.61 / 512;

    % helper: read all distances (in µm) from a sample directory
    getSampleDists = @(sname) readAllDistances( ...
        fullfile(analysisRoot, sname), pixToUm);

    % ---------------------------------------------------------------------
    % Define the series for which we want ΔD
    % ---------------------------------------------------------------------

    % EndoH series at 37C: EndoH−, 1:20, 1:2
    endoSeries(1).label = 'EndoH− (37C, R8L+EB)';
    endoSeries(1).minus = {'R8L 37C EndoH- SNA-', 'EB'};
    endoSeries(1).plus  = {'R8L 37C EndoH- SNA+', 'EB SNA'};

    endoSeries(2).label = 'EndoH 1:20 (E10)';
    endoSeries(2).minus = {'E10'};
    endoSeries(2).plus  = {'E10 SNA'};

    endoSeries(3).label = 'EndoH 1:2 (37C, R8L+E2)';
    endoSeries(3).minus = {'R8L 37C EndoH+ SNA-', 'E2'};
    endoSeries(3).plus  = {'R8L 37C EndoH+ SNA+', 'E2 SNA'};

    % NA series: NB, N10, N2
    naSeries(1).label = 'NB';
    naSeries(1).minus = {'NB'};
    naSeries(1).plus  = {'NB SNA'};

    naSeries(2).label = 'N10';
    naSeries(2).minus = {'N10'};
    naSeries(2).plus  = {'N10 SNA'};

    naSeries(3).label = 'N2';
    naSeries(3).minus = {'N2'};
    naSeries(3).plus  = {'N2 SNA'};

    % ---------------------------------------------------------------------
    % Compute ΔD arrays
    % ---------------------------------------------------------------------
    deltaEndo = computeDeltaForSeries(endoSeries, getSampleDists);
    deltaNA   = computeDeltaForSeries(naSeries,   getSampleDists);

    % ---------------------------------------------------------------------
    % Plot
    % ---------------------------------------------------------------------
    x = 1:3;
    xLabels = {'No enzyme','Low enzyme','High enzyme'};

    figure('Color','w','Name','ΔD trend: EndoH vs NA');
    hold on;

    plot(x, deltaEndo, '-o', ...
        'LineWidth', 1.8, ...
        'MarkerSize', 8);
    plot(x, deltaNA, '-s', ...
        'LineWidth', 1.8, ...
        'MarkerSize', 8);

    xlim([0.8 3.2]);
    xticks(x);
    xticklabels(xLabels);

    ylabel('\DeltaD (mean SNA^+ - mean SNA^-) [\mum]');
    grid on;
    legend({'EndoH series','NA series'}, 'Location','best');
    title('\DeltaD vs enzyme level');

    hold off;
end

% =====================================================================
function deltaVals = computeDeltaForSeries(seriesDefs, getSampleDists)
%COMPUTEDELTAFORSERIES  mean(SNA+) - mean(SNA−) per enzyme level

    n = numel(seriesDefs);
    deltaVals = nan(1,n);

    for i = 1:n
        minusD = [];
        plusD  = [];

        % pool all "minus" samples
        for j = 1:numel(seriesDefs(i).minus)
            minusD = [minusD; getSampleDists(seriesDefs(i).minus{j})]; %#ok<AGROW>
        end

        % pool all "plus" samples
        for j = 1:numel(seriesDefs(i).plus)
            plusD = [plusD; getSampleDists(seriesDefs(i).plus{j})]; %#ok<AGROW>
        end

        if ~isempty(minusD) && ~isempty(plusD)
            muMinus = mean(minusD, 'omitnan');
            muPlus  = mean(plusD,  'omitnan');
            deltaVals(i) = muPlus - muMinus;
        else
            deltaVals(i) = NaN;
        end
    end
end

% =====================================================================
% Re-use the same readAllDistances helper used in your other plotting code
% (if you already have this in another file, you can remove this copy).
function dUm = readAllDistances(sampleDir, pixToUm)
    dUm = [];
    if ~isfolder(sampleDir), return; end
    files = dir(fullfile(sampleDir, '*_distances.xlsx'));
    for iF = 1:numel(files)
        fpath = fullfile(sampleDir, files(iF).name);
        try
            T = readtable(fpath);
        catch
            continue;
        end
        if ~ismember('DistancePixels', T.Properties.VariableNames)
            continue;
        end
        d = T.DistancePixels(:);
        d = d(~isnan(d) & d > 0);
        dUm = [dUm; d * pixToUm]; %#ok<AGROW>
    end
end
