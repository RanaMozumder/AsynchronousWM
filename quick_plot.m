function l = quick_plot(t, y, clr, varargin)
fill([t', fliplr(t')], [mean(y, 1, 'omitnan') + std(y, 0, 1, 'omitnan')./sqrt(size(y, 1)), fliplr(mean(y, 1, 'omitnan') - std(y, 0, 1, 'omitnan')./sqrt(sum(~isnan(y), 1)))], clr, 'FaceAlpha', 0.1, 'EdgeColor','none')
l = plot(t, mean(y, 1, 'omitnan'), 'color', clr, varargin{:});
end