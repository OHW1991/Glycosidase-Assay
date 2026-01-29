function generateFoVFigures(FoVs, outputRoot)
%GENERATEFOVFIGURES  Create 2x4 per-channel+composite figures for each FoV.
%
%   generateFoVFigures(FoVs, 'FoVs_20251130');
%
%   FoVs      : struct array (fields: c1, c2, c3, masterName, conditionName, fovName).
%               c1/c2/c3 can be either full paths (char/string) or image
%               matrices already loaded.
%   outputRoot: top-level folder where figures will be saved
%               (default: 'FoVs_20251130').
%
%   This ADAPTED VERSION:
%     • Only processes FoVs from the date 20251130 when the channel fields
%       are path strings containing "20251130".
%     • If c1/c2/c3 are already matrices (no path info), it will process
%       them regardless of date (no way to detect date).
%
%   For each qualifying FoV, this creates a 2x4 figure:
%
%       RAW / quantitative (absolute 0..65535 scaling):
%         (1,1) c1 RAW    (RED pseudocolor + colorbar)
%         (1,2) c2 RAW    (GREEN pseudocolor + colorbar)
%         (2,1) c3 RAW    (BLUE pseudocolor + colorbar)
%         (2,2) Composite RAW (c1->R, c2->G, c3->B)
%
%       OPTIMIZED / visual:
%         (1,3) c1 OPT    (RED, per-FoV contrast + colorbar)
%         (1,4) c2 OPT    (GREEN, per-FoV contrast EXCEPT SNA- + colorbar)
%         (2,3) c3 OPT    (BLUE, per-FoV contrast + colorbar)
%         (2,4) Composite OPT
%
%   - RAW panels preserve quantitative comparison (0..65535 mapping).
%   - OPT panels enhance visibility.
%   - For SNA- conditions (conditionName contains 'SNA-'), c2 is NOT
%     rescaled in the OPT panels to avoid boosting pure noise.
%
%   SCALE BAR:
%     - Scale bars (5 µm) are drawn on all single-channel panels.
%     - Pixel size is set from 101.61 µm per 512 pixels.

    % ----- default output folder: dedicated for this date -----
    if nargin < 2 || isempty(outputRoot)
        outputRoot = 'FoVs_20251130';
    end

    if ~isfolder(outputRoot)
        mkdir(outputRoot);
    end

    % ==== Microscope scaling (from 101.61 µm / 512 px) ====
    pixelSizeUm  = 101.61 / 512;   % µm per pixel
    barLengthUm  = 5;              % length of scale bar in µm
    % ======================================================

    % Colormaps for channels
    N = 256;
    cmapRed   = [linspace(0,1,N)', zeros(N,1), zeros(N,1)];
    cmapGreen = [zeros(N,1), linspace(0,1,N)', zeros(N,1)];
    cmapBlue  = [zeros(N,1), zeros(N,1), linspace(0,1,N)'];

    nFoVs = numel(FoVs);
    fprintf('Generating FoV figures for %d FoVs (filtering to 20251130 when possible)...\n', nFoVs);

    nUsed = 0;

    for i = 1:nFoVs
        F = FoVs(i);

        % -----------------------------------------------------------------
        % Extract channel "sources" (could be paths or matrices)
        % -----------------------------------------------------------------
        c1src = [];
        c2src = [];
        c3src = [];

        if isfield(F,'c1'), c1src = F.c1; end
        if isfield(F,'c2'), c2src = F.c2; end
        if isfield(F,'c3'), c3src = F.c3; end

        % Basic presence check
        if isempty(c1src) || isempty(c2src) || isempty(c3src)
            warning('Skipping FoV %d (missing c1/c2/c3).', i);
            continue;
        end

        % -----------------------------------------------------------------
        % DATE FILTER: only use FoVs that come from 20251130
        % We can only test this if the sources are path strings.
        % If all channels are matrices (no path info), we process them
        % (cannot tell the date).
        % -----------------------------------------------------------------
        inThisDate = false;
        anyPath    = false;

        for src = {c1src, c2src, c3src}
            s = src{1};
            if ischar(s) || isstring(s)
                anyPath = true;
                if contains(char(s), '20251130')
                    inThisDate = true;
                    break;
                end
            end
        end

        if anyPath && ~inThisDate
            % We have path-based channels but none mentions 20251130 → skip
            continue;
        end

        % -----------------------------------------------------------------
        % Load images (if already matrices, just pass them through)
        % -----------------------------------------------------------------
        try
            I1 = loadChannelImage(c1src);
            I2 = loadChannelImage(c2src);
            I3 = loadChannelImage(c3src);
        catch ME
            warning('Error loading images for FoV %d: %s', i, ME.message);
            continue;
        end

        % Ensure same XY size
        if ~isequal(size(I1,1), size(I2,1), size(I3,1)) || ...
           ~isequal(size(I1,2), size(I2,2), size(I3,2))
            warning('Mismatched image sizes in FoV %d (master=%s, cond=%s, fov=%s). Skipping.', ...
                    i, safeStr(F,'masterName'), safeStr(F,'conditionName'), safeStr(F,'fovName'));
            continue;
        end

        % Detect SNA- from condition name
        condName = safeStr(F,'conditionName');
        isSNAminus = contains(condName, 'SNA-');

        % ---- RAW: scale to [0,1] assuming 16-bit images (0..65535) ----
        I1_raw = castTo01(I1);
        I2_raw = castTo01(I2);
        I3_raw = castTo01(I3);

        composite_raw = cat(3, I1_raw, I2_raw, I3_raw);

        % ---- OPTIMIZED: per-FoV contrast ----
        I1_opt = autoScale01(I1);
        I3_opt = autoScale01(I3);

        % c2 optimization: SNA- keeps RAW scaling
        if isSNAminus
            I2_opt = I2_raw;
        else
            I2_opt = autoScale01(I2);
        end

        composite_opt = cat(3, I1_opt, I2_opt, I3_opt);

        % -----------------------------------------------------------------
        % Prepare output path: FoVs_20251130/master/condition
        % -----------------------------------------------------------------
        masterName = safeStr(F,'masterName');
        fovName    = safeStr(F,'fovName');

        outDir = fullfile(outputRoot, masterName, condName);
        if ~isfolder(outDir)
            mkdir(outDir);
        end

        % Safe FoV name for file
        safeFovName = regexprep(fovName, '[^\w\-]', '_');
        if isempty(safeFovName)
            safeFovName = sprintf('FoV_%d', i);
        end
        outPDF = fullfile(outDir, [safeFovName '.pdf']);
        outJPG = fullfile(outDir, [safeFovName '.jpg']);

        % Create figure off-screen (wider for 2x4)
        hFig = figure('Visible','off', ...
                      'Position',[100 100 1600 800], ...
                      'Color','w');

        % ===================== ROW 1 =====================

        % (1,1) c1 RAW
        ax = subplot(2,4,1);
        imagesc(ax, I1_raw);
        axis(ax, 'image'); axis(ax, 'off');
        colormap(ax, cmapRed);
        caxis(ax, [0 1]);
        title(ax, sprintf('c1 RAW (RED)\n%s | %s | %s', ...
            masterName, condName, fovName), ...
            'Interpreter','none','FontSize',8);
        addColorbarRaw(ax);
        addScaleBar(ax, pixelSizeUm, barLengthUm, [1 1 1]);

        % (1,2) c2 RAW
        ax = subplot(2,4,2);
        imagesc(ax, I2_raw);
        axis(ax, 'image'); axis(ax, 'off');
        colormap(ax, cmapGreen);
        caxis(ax, [0 1]);
        title(ax, 'c2 RAW (GREEN)', 'FontSize',8);
        addColorbarRaw(ax);
        addScaleBar(ax, pixelSizeUm, barLengthUm, [1 1 1]);

        % (1,3) c1 OPT
        ax = subplot(2,4,3);
        imagesc(ax, I1_opt);
        axis(ax, 'image'); axis(ax, 'off');
        colormap(ax, cmapRed);
        caxis(ax, [0 1]);
        title(ax, 'c1 OPT (RED)', 'FontSize',8);
        addColorbarOpt(ax);
        addScaleBar(ax, pixelSizeUm, barLengthUm, [1 1 1]);

        % (1,4) c2 OPT
        ax = subplot(2,4,4);
        imagesc(ax, I2_opt);
        axis(ax, 'image'); axis(ax, 'off');
        colormap(ax, cmapGreen);
        caxis(ax, [0 1]);
        if isSNAminus
            title(ax, 'c2 OPT (GREEN, SNA-: no rescale)', 'FontSize',8);
        else
            title(ax, 'c2 OPT (GREEN, rescaled)', 'FontSize',8);
        end
        addColorbarOpt(ax);
        addScaleBar(ax, pixelSizeUm, barLengthUm, [1 1 1]);

        % ===================== ROW 2 =====================

        % (2,1) c3 RAW
        ax = subplot(2,4,5);
        imagesc(ax, I3_raw);
        axis(ax, 'image'); axis(ax, 'off');
        colormap(ax, cmapBlue);
        caxis(ax, [0 1]);
        title(ax, 'c3 RAW (BLUE)', 'FontSize',8);
        addColorbarRaw(ax);
        addScaleBar(ax, pixelSizeUm, barLengthUm, [1 1 1]);

        % (2,2) Composite RAW
        ax = subplot(2,4,6);
        imshow(composite_raw, 'Parent', ax);
        axis(ax, 'image'); axis(ax, 'off');
        title(ax, 'Composite RAW (c1→R, c2→G, c3→B)', 'FontSize',8);

        % (2,3) c3 OPT
        ax = subplot(2,4,7);
        imagesc(ax, I3_opt);
        axis(ax, 'image'); axis(ax, 'off');
        colormap(ax, cmapBlue);
        caxis(ax, [0 1]);
        title(ax, 'c3 OPT (BLUE)', 'FontSize',8);
        addColorbarOpt(ax);
        addScaleBar(ax, pixelSizeUm, barLengthUm, [1 1 1]);

        % (2,4) Composite OPT
        ax = subplot(2,4,8);
        imshow(composite_opt, 'Parent', ax);
        axis(ax, 'image'); axis(ax, 'off');
        if isSNAminus
            title(ax, 'Composite OPT (G from RAW, SNA-)', 'FontSize',8);
        else
            title(ax, 'Composite OPT (all channels rescaled)', 'FontSize',8);
        end

        % ---- Save ----
        try
            exportgraphics(hFig, outPDF, 'ContentType','vector');
        catch ME
            warning('Failed to save PDF for %s: %s', safeFovName, ME.message);
        end

        try
            exportgraphics(hFig, outJPG, 'Resolution',300);
        catch ME
            warning('Failed to save JPG for %s: %s', safeFovName, ME.message);
        end

        close(hFig);

        nUsed = nUsed + 1;
        fprintf('Saved FoV (20251130 candidate): %s (PDF+JPG)\n', safeFovName);
    end

    fprintf('Done generating FoV figures. Used %d FoVs.\n', nUsed);
end

% =====================================================================
% Helper: load channel image from path or pass-through matrix
% =====================================================================
function I = loadChannelImage(src)
    if ischar(src) || isstring(src)
        I = imread(src);
    elseif isnumeric(src)
        I = src;
    else
        error('Unsupported channel source type: %s', class(src));
    end
end

% =====================================================================
% Helper: safe string extraction from struct field
% =====================================================================
function s = safeStr(S, fieldName)
    if isfield(S, fieldName) && ~isempty(S.(fieldName))
        val = S.(fieldName);
        if isstring(val) || ischar(val)
            s = char(val);
        else
            s = sprintf('%s', mat2str(val));
        end
    else
        s = '';
    end
end

% -------------------------------------------------------------------------
function J = castTo01(I)
%CASTTO01  Map integer images to [0,1] with a 16-bit assumption when possible.

    if isa(I, 'uint16')
        J = double(I) / 65535;
    elseif isinteger(I)
        J = double(I) / double(intmax(class(I)));
    else
        J = mat2gray(I);
    end
end

% -------------------------------------------------------------------------
function J = autoScale01(I)
%AUTOSCALE01  Per-FoV contrast stretching to [0,1].

    I = double(I);
    minVal = min(I(:));
    maxVal = max(I(:));

    if maxVal <= minVal
        J = zeros(size(I));
    else
        J = (I - minVal) / (maxVal - minVal);
    end
end

% -------------------------------------------------------------------------
function addScaleBar(ax, pixelSizeUm, barLengthUm, color)
%ADDSCALEBAR  Draw a scale bar on an image axis without changing layout.

    if pixelSizeUm <= 0 || barLengthUm <= 0
        return;
    end

    axes(ax); %#ok<LAXES>
    hold(ax, 'on');

    xL = get(ax, 'XLim');
    yL = get(ax, 'YLim');

    widthPx  = diff(xL);
    heightPx = diff(yL);

    barPx = barLengthUm / pixelSizeUm;
    if barPx > 0.5*widthPx
        barPx = 0.5 * widthPx;
    end

    marginX = 0.05 * widthPx;
    marginY = 0.05 * heightPx;

    xStart = xL(1) + marginX;
    xEnd   = xStart + barPx;
    yBar   = yL(2) - marginY;

    plot(ax, [xStart, xEnd], [yBar, yBar], 'Color', color, 'LineWidth', 2);

    labelStr = sprintf('%g \\mum', barLengthUm);
    text(ax, (xStart + xEnd)/2, yBar - 0.02*heightPx, labelStr, ...
        'Color', color, 'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', 'FontSize',8);

    hold(ax, 'off');
end

% -------------------------------------------------------------------------
function addColorbarRaw(ax)
%ADDCOLORBARRAW  Add a 16-bit intensity colorbar without shrinking the axis.

    origPos = get(ax, 'Position');
    cb = colorbar(ax);
    cb.Location = 'eastoutside';
    tickVals   = [0 0.25 0.5 0.75 1];
    tickLabels = {'0','16384','32768','49152','65535'};
    set(cb, 'Ticks', tickVals, 'TickLabels', tickLabels);
    ylabel(cb, 'Intensity (16-bit)');
    set(ax, 'Position', origPos);
end

% -------------------------------------------------------------------------
function addColorbarOpt(ax)
%ADDCOLORBAROPT  Add a normalized intensity colorbar without shrinking axis.

    origPos = get(ax, 'Position');
    cb = colorbar(ax);
    cb.Location = 'eastoutside';
    set(cb, 'Ticks', [0 0.5 1], 'TickLabels', {'0','0.5','1'});
    ylabel(cb, 'Norm. intensity (0–1)');
    set(ax, 'Position', origPos);
end
