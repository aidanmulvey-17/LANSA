function h = plotShaded(x, y, err, color, varargin)
%PLOTSHADED Line with a shaded error band (y - err to y + err).
%   h = plotShaded(x, y, err, color, ...) draws the band (alpha 0.25) and
%   the line; extra arguments are passed to plot (e.g. 'DisplayName').
x = x(:)';  y = y(:)';  err = err(:)';
fill([x, fliplr(x)], [y + err, fliplr(y - err)], color, 'EdgeColor', 'none', ...
    'FaceAlpha', 0.25, 'HandleVisibility', 'off');
hold on;
h = plot(x, y, 'Color', color, 'LineWidth', 2, varargin{:});
end
