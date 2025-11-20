clc; clearvars; %close all;
%% session and parameters
sessions = readtable("D:\NeuropixelData\CellReportsDataCode\Data\NeuropixelData\SessionInfo.xlsx", 'Sheet', 'EligibleSessions');
least_n = 5;
statecode_threshold = 7;%[5, 6]; 
t_range = [0.5 3.5]; 
pop = 0;
%% create spike populations
for ses=1:height(sessions)
    ses_name = cell2mat(sessions{ses, 1}) 
    load(['D:\NeuropixelData\CellReportsDataCode\Data\NeuropixelData\NeuronFiles\', ses_name, '_1_sparse.mat'])
    cue_neurons = readtable(['D:\NeuropixelData\CellReportsDataCode\Data\NeuropixelData\SelectiveNeurons\', ses_name,'.xlsx'], ...
        'Sheet', 'cuedelay_selective');
    for cl=1:8
        current = cue_neurons(cue_neurons{:, 9} == cl, :);
        if height(current) < least_n
            continue
        else
            % for correct trials
            target_trials = find([MatData.trials.Statecode] == statecode_threshold & [MatData.trials.Class] ==cl); 
%             target_trials = find([MatData.trials.Statecode] == statecode_threshold); 
            % for error trials
%             target_trials = find(([MatData.trials.Statecode] == 5 | [MatData.trials.Statecode] == 6) & [MatData.trials.Class] ==cl); 
            if length(target_trials)>4
                pop = pop +1;
            else
                continue
            end
        end
        [~, ind] = ismember(cellfun(@str2double, current.Cluster)', MatData.cids);
        corr_tr = 0;
        spike_times = cell(numel(ind), numel(target_trials));
        for tr = target_trials
            corr_tr = corr_tr + 1;
            t_shift = double(MatData.trials(tr).photodiode_on_event)/MatData.sample_rate;
            t_delay = t_range ;
            ss_mat = MatData.trials(tr).ss(:, ind);
            for k = 1:size(ss_mat,2)
                st = (find(ss_mat(:, k))/MatData.sample_rate)' - t_shift;
                st_delay = st(st>t_delay(1) & st<=t_delay(2));
                spike_times{k, corr_tr}=st_delay;
            end
            clear ss_mat
        end
        spk(pop) = {spike_times};
    end
    clear MatData   
end
clearvars -except spk
%% Step 2: Raster Plot Function
% figure;
duration = 3.5;
% plot_raster(spk{1}, duration, 'Without Silent Periods');

%% Calculate ISI distributions
num_shuffles = 10000;
if isempty(gcp('nocreate'))
    parpool; % Adjust the number of workers as needed
end
for pop = 1:length(spk)
    fprintf('%d ',pop)
    spk_pop   = spk{pop};
    isi_pop   = [];
    pop_trials = size(spk_pop, 2);
    pop_neurons = size(spk_pop, 1);
    for trial = 1:pop_trials
        spk_temp   = [];
        for neuron = 1:pop_neurons
            spk_temp   = [spk_temp, spk_pop{neuron, trial}];
        end
        spk_temp    = sort(spk_temp);
        freq(trial) = numel(spk_temp)/duration;
        isi_pop     = [isi_pop, diff(spk_temp)];

    end
    max_isi_emp(pop) = max(isi_pop);
    isi_emp(pop)     = {isi_pop};
    freq_emp(pop)    = mean(freq);

    %null
    max_isi_pop_null    = zeros(num_shuffles, 1);
    freq_null_shuffle   = zeros(num_shuffles, 1);

    parfor s = 1:num_shuffles
        isi_pop_null_local   = [];
        freq_local   = zeros(pop_trials, 1);

        for trial = 1:pop_trials
            if pop_trials>=pop_neurons
                perm_trials = randperm(pop_trials, pop_neurons);
            else
                d = floor(pop_neurons/pop_trials);
                r = rem(pop_neurons, pop_trials);
                perm_trials = [];
                for i = 1:d
                    perm_trials = [perm_trials, randperm(pop_trials, pop_trials)];
                end
                perm_trials = [perm_trials, randperm(pop_trials, r)];
            end
            perm_trials = perm_trials(randperm(numel(perm_trials)));
            spk_temp   = [];
            for neuron = 1:pop_neurons
                spk_temp   = [spk_temp, spk_pop{neuron, perm_trials(neuron)}];
            end
            spk_temp           = sort(spk_temp);
            freq_local(trial)  = numel(spk_temp) / duration;
            isi_pop_null_local = [isi_pop_null_local, diff(spk_temp)];
        end
        if ~isempty(isi_pop_null_local)
            max_isi_pop_null(s)  = max(isi_pop_null_local);
        else
            max_isi_pop_null(s) = nan;
        end
        freq_null_shuffle(s) = mean(freq_local);
    end
    isi_null{pop}   = max_isi_pop_null;

    max_isi_null(pop,:)   = max_isi_pop_null';
    freq_null(pop,:)   = freq_null_shuffle';

end
%% fit linear model
slope_null = nan(num_shuffles, 1);
intercept_null = nan(num_shuffles, 1);

for s=1:num_shuffles
    c = polyfit(log10(freq_null(:,s)), log10(max_isi_null(:,s)), 1);

    slope_null(s)     = c(1);
    intercept_null(s) = c(2);
    clear c
end
c = polyfit(log10(freq_emp), log10(max_isi_emp), 1);

slope_emp     = c(1);
intercept_emp = c(2);

%% 
fig = figure('Units','inches', 'Position',[0 0 5 2.5]);  % Total size ~5x2.5 inches
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% --------- Slope Histogram ---------
nexttile;
h1 = histogram(slope_null, 'Normalization', 'probability', ...
    'FaceColor', [0.2 0.4 1], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
hold on;
xline(slope_emp, 'k', 'LineWidth', 2);
legend('Shuffled slopes', 'Empirical slope', 'Location', 'best');
xlabel('Slope'); ylabel('Probability');
set(gca, 'FontSize', 12)

% --------- Intercept Histogram ---------
nexttile;
h2 = histogram(intercept_null, 'Normalization', 'probability', ...
    'FaceColor', [0.2 0.4 1], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
hold on;
xline(intercept_emp, 'k', 'LineWidth', 2);
legend('Shuffled intercepts', 'Empirical intercept', 'Location', 'best');
xlabel('Intercept'); ylabel('Probability');
set(gca, 'FontSize', 12)


%% 
figure('Units','inches', 'Position',[0 0 4 4]); hold on
plot(log10(freq_emp), log10(max_isi_emp), 'o', 'MarkerSize',5, 'MarkerEdgeColor','k')
plot(0:0.01:3, polyval(c, 0:0.01:3), 'k-')
plot_lm_ci([slope_null, intercept_null], 0:0.01:3, [0.025 0.975])
% xlim([0.8 2.6])
xlabel('log_1_0(frequency)')
ylabel('log_1_0(max-isi)')
legend('Empirical data points', 'Empirical fit', 'Shuffled')
set(gca, 'FontSize', 12)

%%


