function p_value = perm_test_single_psth(baseline_rate, task_rate)

% Compute observed mean firing rate in task epoch
observed_mean_task = mean(task_rate);

% Permutation test
all_bins = [baseline_rate, task_rate];
labels = [zeros(1, length(baseline_rate)), ones(1, length(task_rate))];
n_iter = 10000;
shuffled_means = zeros(n_iter, 1);

for i = 1:n_iter
    shuffled_labels = labels(randperm(length(labels)));
    shuffled_task = all_bins(shuffled_labels == 1);
    shuffled_means(i) = mean(shuffled_task);
end

% Compute p-value
p_value = sum(shuffled_means >= observed_mean_task) / n_iter;
end


