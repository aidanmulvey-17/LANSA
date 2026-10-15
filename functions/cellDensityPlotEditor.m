function modifiedData = cellDensityPlotEditor(initialData, szData, depth, layerBounds, cellTypes)
%CELLDENSITYPLOTEDITOR Interactive editor for laminar cell-density profiles.
%   modifiedData = cellDensityPlotEditor(initialData, szData, depth, layerBounds, cellTypes)
%   opens a window in which the density profile of each cell type can be
%   reshaped by clicking a point and dragging it left or right.
%
%   initialData : bins x 4 control densities (grey solid reference line)
%   szData      : bins x 4 schizophrenia-modeled densities (grey dashed line)
%   depth       : bins x 1 laminar position of each bin
%   layerBounds : laminar boundaries from pia to white matter (layer lines)
%   cellTypes   : names of the four columns (default {'PV','CB','CR','NG'})
%
%   Controls
%     Cell type dropdown   choose which profile to edit
%     Schizophrenia box    set every cell type to the schizophrenia-reported
%                          densities (untick to return to control densities);
%                          points can still be dragged afterwards
%     Reset                discard all edits (control densities)
%     Submit               return the edited bins x 4 matrix (closing the
%                          window also returns the current values)

if nargin < 5 || isempty(cellTypes), cellTypes = {'PV', 'CB', 'CR', 'NG'}; end
depth = depth(:);
dataMatrix  = initialData;
currentType = 1;
draggingBin = 0;

fig = uifigure('Name', 'Interactive Cortical Density Editor', ...
    'Position', [100 100 600 740], ...
    'WindowButtonMotionFcn', @dragPoint, ...
    'WindowButtonUpFcn', @stopDrag, ...
    'CloseRequestFcn', @(src, ~) uiresume(src));

% ---- Controls ------------------------------------------------------------
uilabel(fig, 'Position', [50 690 100 22], 'Text', 'Select Cell Type:');
uidropdown(fig, 'Items', cellTypes, 'Position', [150 690 150 22], ...
    'ValueChangedFcn', @changeCellType);
uibutton(fig, 'Text', 'Submit Final Values & Plot', 'Position', [350 686 200 30], ...
    'ButtonPushedFcn', @(~, ~) uiresume(fig));
szBox = uicheckbox(fig, 'Text', 'Use schizophrenia-reported densities (all cell types)', ...
    'Position', [50 655 330 22], 'ValueChangedFcn', @toggleSchizophrenia);
uibutton(fig, 'Text', 'Reset', 'Position', [450 652 100 26], ...
    'ButtonPushedFcn', @resetData);

% ---- Axes ----------------------------------------------------------------
ax = uiaxes(fig, 'Position', [50 50 500 590]);
disableDefaultInteractivity(ax);   % so that click-and-drag edits the points
xMax = 1.25 * max([initialData(:); szData(:)]);
ax.YDir  = 'reverse';              % pia at the top
ax.YLim  = [depth(1) - 0.05, depth(end) + 0.05];
ax.YTick = depth(1:2:end);
ax.XLim  = [0 xMax];
grid(ax, 'on');
ylabel(ax, 'Cortical depth (relative to mid-L4)');
xlabel(ax, 'Cell density');
title(ax, 'Click and drag points left or right to adjust density');
hold(ax, 'on');

innerBounds = layerBounds(2:end-1);
layerMids   = (layerBounds(1:end-1) + layerBounds(2:end)) / 2;
yline(ax, innerBounds, '--', 'HandleVisibility', 'off');
text(ax, 0.9 * xMax * ones(size(layerMids)), layerMids, ...
    {'L1', 'L2', 'L3', 'L4', 'L5', 'L6'}, 'FontWeight', 'bold');

% Reference lines (drawn first so they stay in the background)
hRefPlot = plot(ax, initialData(:, currentType), depth, '-', ...
    'Color', [0.6 0.6 0.6], 'LineWidth', 1.5);
hSzPlot  = plot(ax, szData(:, currentType), depth, '--', ...
    'Color', [0.6 0.6 0.6], 'LineWidth', 1.5);

% Editable line; only this line starts a drag
hPlot = plot(ax, dataMatrix(:, currentType), depth, '-o', 'MarkerSize', 8, ...
    'MarkerFaceColor', [0 0.4470 0.7410], 'LineWidth', 2);
hPlot.ButtonDownFcn = @startDrag;

% ---- Callbacks -------------------------------------------------------------
    function startDrag(~, event)
        % Select the bin closest to the clicked point
        [~, draggingBin] = min(abs(depth - event.IntersectionPoint(2)));
    end

    function dragPoint(~, ~)
        if draggingBin > 0
            pt = ax.CurrentPoint;
            dataMatrix(draggingBin, currentType) = max(0, min(ax.XLim(2), pt(1, 1)));
            hPlot.XData = dataMatrix(:, currentType);
        end
    end

    function stopDrag(~, ~)
        draggingBin = 0;
    end

    function changeCellType(src, ~)
        currentType = find(strcmp(src.Items, src.Value));
        hPlot.XData    = dataMatrix(:, currentType);
        hRefPlot.XData = initialData(:, currentType);
        hSzPlot.XData  = szData(:, currentType);
    end

    function toggleSchizophrenia(src, ~)
        if src.Value
            dataMatrix = szData;
        else
            dataMatrix = initialData;
        end
        hPlot.XData = dataMatrix(:, currentType);
    end

    function resetData(~, ~)
        dataMatrix = initialData;
        szBox.Value = false;
        hPlot.XData = dataMatrix(:, currentType);
    end

% ---- Wait for Submit (or window close) -------------------------------------
uiwait(fig);
modifiedData = dataMatrix;
if isvalid(fig)
    delete(fig);
end
end
