clear all; clearvars%close all; clc;

% Parameters
num_neurons = 2:20;       % Number of neurons
num_trials = 10:20;        % Number of trials
duration = 3000;        % Duration of each trial in ms
silent_period_dur = 100; % Duration of each silent period (ms)
num_pop = 50;           % Number of populations to simulate
num_shuffles = 10000;   % Number of shuffles to generate null
min_silents = 1; max_silents = 3; % Range for number of silent periods per trial
%% Generate surrogate data
for pop = 1:num_pop
    pop_neuron = num_neurons(randi(numel(num_neurons)));
    pop_trials = num_trials(randi(numel(num_trials)));

    spike_times_no_silence = cell(pop_neuron, pop_trials);
    spike_times_with_silence = cell(pop_neuron, pop_trials);

    % Random number of silent periods per trial
    num_silent_periods = randi([min_silents, max_silents], [1, pop_trials]);
    silent_windows = cell(1, pop_trials);

    for trial = 1:pop_trials
        silent_onsets = sort(randperm(duration - silent_period_dur, num_silent_periods(trial)));
        silent_windows{trial} = [silent_onsets' silent_onsets' + silent_period_dur];
    end

    for neuron = 1:pop_neuron
        baseline_FR = randi(20);  % Random FR per neuron

        for trial = 1:pop_trials
            % Generate spike train with no silence
            spike_times = sort(randpoisson(baseline_FR, duration));
            spike_times_no_silence{neuron, trial} = spike_times;

            % Generate spike train with silence
            spike_times = sort(randpoisson(baseline_FR, duration));
            % Remove spikes in all silent windows
            for w = 1:size(silent_windows{trial}, 1)
                win = silent_windows{trial}(w, :);
                spike_times = spike_times(spike_times < win(1) | spike_times > win(2));
            end
            spike_times_with_silence{neuron, trial} = spike_times;
        end
    end

    spk(pop) = {spike_times_no_silence};
    spk_silence(pop) = {spike_times_with_silence};

    clear spike_times_with_silence spike_times_no_silence
end

%% Step 2: Raster Plot Function
figure('Units','inches', 'Position',[0 0 8 2.5]);
plot_raster(spk{1}, duration, 'Without Silent Periods'); axis tight
xticklabels([0.5 1 1.5 2 2.5 3])
figure('Units','inches', 'Position',[0 0 8 2.5]);
plot_raster(spk_silence{1}, duration, 'With Silent Periods'); axis tight
xticklabels([0.5 1 1.5 2 2.5 3])

%% Calculate ISI distributions
if isempty(gcp('nocreate'))
    parpool; % Adjust the number of workers as needed
end
for pop = 1:num_pop
    fprintf('%d ',pop)
    spk_pop   = spk{pop};
    spk_pop_s = spk_silence{pop};
    isi_pop   = [];
    isi_pop_s = [];
    pop_trials = size(spk_pop, 2);
    pop_neurons = size(spk_pop, 1);
    for trial = 1:pop_trials
        spk_temp   = [];
        spk_temp_s = [];
        for neuron = 1:pop_neurons
            spk_temp   = [spk_temp, spk_pop{neuron, trial}];
            spk_temp_s = [spk_temp_s, spk_pop_s{neuron, trial}];
        end
        spk_temp   = sort(spk_temp);
        spk_temp_s = sort(spk_temp_s);

        freq(trial)   = numel(spk_temp)/(duration/1000);
        freq_s(trial) = numel(spk_temp_s)/(duration/1000);

        isi_pop   = [isi_pop diff(spk_temp)];
        isi_pop_s = [isi_pop_s diff(spk_temp_s)];

    end
    max_isi_emp(pop)   = max(isi_pop);
    max_isi_emp_s(pop) = max(isi_pop_s);

    isi_emp(pop)   = {isi_pop};
    isi_emp_s(pop) = {isi_pop_s};

    freq_emp(pop)   = mean(freq);
    freq_emp_s(pop) = mean(freq_s);

    clear freq freq_s
    %null
    max_isi_pop_null    = zeros(num_shuffles, 1);
    max_isi_pop_null_s  = zeros(num_shuffles, 1);
    freq_null_shuffle   = zeros(num_shuffles, 1);
    freq_null_shuffle_s = zeros(num_shuffles, 1);

    parfor s = 1:num_shuffles
        isi_pop_null_local   = [];
        isi_pop_null_s_local = [];
        freq_local   = zeros(pop_trials, 1);
        freq_s_local = zeros(pop_trials, 1);

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
            spk_temp_s = [];

            for neuron = 1:pop_neurons
                spk_temp   = [spk_temp, spk_pop{neuron, perm_trials(neuron)}];
                spk_temp_s = [spk_temp_s, spk_pop_s{neuron, perm_trials(neuron)}];
            end

            spk_temp   = sort(spk_temp);
            spk_temp_s = sort(spk_temp_s);

            freq_local(trial)   = numel(spk_temp) / (duration/1000);
            freq_s_local(trial) = numel(spk_temp_s) / (duration/1000);

            isi_pop_null_local   = [isi_pop_null_local, diff(spk_temp)];
            isi_pop_null_s_local = [isi_pop_null_s_local, diff(spk_temp_s)];
        end

        max_isi_pop_null(s)   = max(isi_pop_null_local);
        max_isi_pop_null_s(s) = max(isi_pop_null_s_local);

        freq_null_shuffle(s)   = mean(freq_local);
        freq_null_shuffle_s(s) = mean(freq_s_local);
    end

    % After parfor
    isi_null{pop}   = max_isi_pop_null;
    isi_null_s{pop} = max_isi_pop_null_s;

    max_isi_null(pop,:)   = max_isi_pop_null';
    max_isi_null_s(pop,:) = max_isi_pop_null;

    freq_null(pop,:)   = freq_null_shuffle';
    freq_null_s(pop,:) = freq_null_shuffle_s';

end
%% fit linear model
slope_null = nan(num_shuffles, 1);
slope_null_s = nan(num_shuffles, 1);

intercept_null = nan(num_shuffles, 1);
intercept_null_s = nan(num_shuffles, 1);

for s=1:num_shuffles
    c     = polyfit(log10(freq_null(:,s)), log10(max_isi_null(:,s)), 1);
    c_s   = polyfit(log10(freq_null_s(:,s)), log10(max_isi_null_s(:,s)), 1);

    slope_null(s) = c(1);
    slope_null_s(s) = c_s(1);

    intercept_null(s) = c(2);
    intercept_null_s(s) = c_s(2);

    clear c c_s
end
c     = polyfit(log10(freq_emp), log10(max_isi_emp), 1);
c_s   = polyfit(log10(freq_emp_s), log10(max_isi_emp_s), 1);

slope_emp = c(1);
slope_emp_s = c_s(1);

intercept_emp = c(2);
intercept_emp_s = c_s(2);

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
sgtitle('Without Silent Periods')

%%
fig = figure('Units','inches', 'Position',[0 0 5 2.5]);  % Total size ~5x2.5 inches
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% --------- Slope Histogram ---------
nexttile;
h1 = histogram(slope_null_s, 'Normalization', 'probability', ...
    'FaceColor', [0.5 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
hold on;
xline(slope_emp_s, 'k', 'LineWidth', 2);
legend('Shuffled slopes', 'Empirical slope', 'Location', 'best');
xlabel('Slope'); ylabel('Probability');
set(gca, 'FontSize', 12)

% --------- Intercept Histogram ---------
nexttile;
h2 = histogram(intercept_null_s, 'Normalization', 'probability', ...
    'FaceColor', [0.5 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
hold on;
xline(intercept_emp_s, 'k', 'LineWidth', 2);
legend('Shuffled intercepts', 'Empirical intercept', 'Location', 'best');
xlabel('Intercept'); ylabel('Probability');
set(gca, 'FontSize', 12)
sgtitle('With Silent Periods')

%% 
figure('Units','inches', 'Position',[0 0 4 4]); hold on
plot(log10(freq_emp), log10(max_isi_emp), 'o', 'MarkerSize',5)
plot(0.8:0.01:2.6, polyval(c, 0.8:0.01:2.6), 'b-')
plot_lm_ci([slope_null, intercept_null], 0.8:0.01:2.6, [0.025 0.975])
xlim([0.8 2.6])
xlabel('log_1_0(frequency)')
ylabel('log_1_0(max-isi)')
legend('Empirical data points', 'Empirical fit', 'Shuffled')
set(gca, 'FontSize', 12)
sgtitle('Without Silent Periods')
%%
figure('Units','inches', 'Position',[0 0 4 4]); hold on
plot(log10(freq_emp_s), log10(max_isi_emp_s), 'o', 'MarkerSize',5, ...
    'MarkerFaceColor', 'none', 'MarkerEdgeColor',[0.5 0 0])
plot(0:0.01:3, polyval(c_s, 0:0.01:3), 'Color', [0.5 0 0], 'LineWidth', 2)
plot_lm_ci([slope_null_s, intercept_null_s], 0:0.01:3, [0.025 0.975])
xlim([0.8 2.6])
xlabel('log_1_0(frequency)')
ylabel('log_1_0(max-isi)')
legend('Empirical data points', 'Empirical fit', 'Shuffled')
set(gca, 'FontSize', 12)
sgtitle('With Silent Periods')


%% Function Definitions

function spike_times = randpoisson(rate, duration)
% Generate Poisson-distributed spike times
num_spikes = poissrnd(rate * (duration / 1000)); % Expected number of spikes
spike_times = sort(randi([1, duration], 1, num_spikes)); % Random spike times
end

function plot_raster(spike_data, duration, title_text)
% Plot raster for spike data
hold on;
num_neurons = size(spike_data, 1);
num_trials = size(spike_data, 2);
for trial = 1:num_trials
    for neuron = 1:num_neurons
        spikes = spike_data{neuron, trial};
        % Add spikes as vertical lines
        for t = 1:length(spikes)
            line([spikes(t), spikes(t)], [trial-0.45, trial+0.45], 'Color', 'k');
        end
    end
end
xlabel('Time (s)');
ylabel('#trials');
xlim([0, duration]);
%     ylim([0, num_trials + 1]);
ylim([0, 20 + 1]);
title(title_text);
set(gca, 'YDir', 'reverse'); % Invert y-axis for typical raster orientation
set(gca, 'Fontsize', 12)
hold off;
end
