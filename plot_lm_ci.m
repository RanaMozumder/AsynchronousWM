function h=plot_lm_ci(lms, xs, ci)
ys = lms * [xs; ones(size(xs))];
mu = mean(ys, 1);
cis = prctile(ys, ci * 100, 1);
fill([xs, fliplr(xs)], [cis(1,:), fliplr(cis(2,:))], ...
    [0 0 1], 'EdgeColor', 'none', 'FaceAlpha',0.3);  % Light gray shade
h=plot(xs, mu, '-b', 'LineWidth', 1);
% plot(xs, mu, '-k');
% plot(xs, cis(1, :), '--k');
% plot(xs, cis(2, :), '--k');
end