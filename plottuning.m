function [gauss_fit_on, gauss_fit_off, rotated_on, rotated_off] = plottuning(data_on, data_off, data_on_null, data_off_null, on, off, oncolor, offcolor)
on_mat_avg = data_on;
on_mat_null_avg = data_on_null;
off_mat_avg = data_off;
off_mat_null_avg = data_off_null;
for n = 1: size(on_mat_avg, 1)
    [~, maxind] = max(sum([on_mat_avg(n,:);off_mat_avg(n,:)], 1));
    rotated_on_mat(n, :) = circular_rotate(on_mat_avg(n,:), maxind);
    rotated_off_mat(n, :) = circular_rotate(off_mat_avg(n,:), maxind);

%     rotated_on_mat_null(n, :) = circular_rotate(on_mat_null_avg(n,:), maxind);
%     rotated_off_mat_null(n, :) = circular_rotate(off_mat_null_avg(n,:), maxind);
end
for n = 1: size(on_mat_avg, 1)
    for shuff=1:size(on_mat_null_avg,1)
        temp_on=squeeze(on_mat_null_avg(shuff, n,:))';
        temp_off=squeeze(off_mat_null_avg(shuff,n,:))';
        [~, maxind] = max(sum([temp_on;temp_off], 1));
        on_shuff(shuff,:)  = circular_rotate(temp_on, maxind);
        off_shuff(shuff,:) = circular_rotate(temp_off, maxind);
    end
    rotated_on_mat_null(n, :)  = mean(on_shuff, 1);
    rotated_off_mat_null(n, :) = mean(off_shuff, 1);
end

% mean_on  = mean(rotated_on_mat, 1, 'omitnan');
% sem_on   = std(rotated_on_mat, 0, 1, 'omitnan') ./ sqrt(size(rotated_on_mat, 1));
% mean_off = mean(rotated_off_mat, 1, 'omitnan');
% sem_off  = std(rotated_off_mat, 0, 1, 'omitnan') ./ sqrt(size(rotated_off_mat, 1));

rotated_on  = rotated_on_mat-rotated_on_mat_null;
rotated_off = rotated_off_mat-rotated_off_mat_null;
mean_on  = mean(rotated_on_mat-rotated_on_mat_null, 1, 'omitnan');
sem_on   = std(rotated_on_mat-rotated_on_mat_null, 0, 1, 'omitnan') ./ sqrt(size(rotated_on_mat-rotated_on_mat_null, 1));
mean_off = mean(rotated_off_mat-rotated_off_mat_null, 1, 'omitnan');
sem_off  = std(rotated_off_mat-rotated_off_mat_null, 0, 1, 'omitnan') ./ sqrt(size(rotated_off_mat-rotated_off_mat_null, 1));

x = -4:4;
% p1 = scatter(x, mean_on, 10, [1, 0.6, 0], 'filled'); hold on;
p1 = scatter(x, mean_on, 20, oncolor, 'filled'); hold on;
p2 = scatter(x, mean_off, 20, offcolor, 'filled'); hold on;

% errorbar(x, mean_on, sem_on, 'Color', [1, 0.6, 0], 'LineStyle', 'none', 'CapSize', 0);
errorbar(x, mean_on, sem_on, 'Color', oncolor, 'LineStyle', 'none', 'CapSize', 0, 'LineWidth',1);
errorbar(x, mean_off, sem_off, 'Color', offcolor, 'LineStyle', 'none', 'CapSize', 0, 'LineWidth',1);
% DOG(mean_on, 'r')
% singleGaussian(mean_on, [1, 0.6, 0])
gauss_fit_on = singleGaussian(mean_on, oncolor);
% DOG(mean_off, 'b')
gauss_fit_off = singleGaussian(mean_off, offcolor);



xlim([-5 5]);
xticks([-4:4])
xticklabels({'-180', '-135', '-90', '-45', '0', '45', '90', '135', '180'})
xlabel('Angular distance from preferred cue (degree)');
ylabel('Normalized mean firing rate');
legend([p1, p2], sprintf('On = %d', on), sprintf('Off = %d', off))
end