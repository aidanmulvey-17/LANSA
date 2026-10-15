function styleAxes(ax)
%STYLEAXES Common axis styling for all LANSA figures.
if nargin < 1, ax = gca; end
set(ax, 'Box', 'off', 'TickDir', 'out', 'Color', 'none', 'LineWidth', 1.2, ...
    'TickLength', [0.02 0.025], 'LooseInset', max(get(ax, 'TightInset'), 0.02));
end
