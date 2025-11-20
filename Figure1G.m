clear all; close all; clc;

sessions = readtable('D:\NeuropixelData\CellReportsDataCode\Data\NeuropixelData\SessionInfo.xlsx', sheet = 'EligibleSessions');

neuron_type = 'single-unit'; % mua or single-unit or delay-unit

for ses = 1%:height(sessions)
    session = cell2mat(sessions{ses, 1})
    load(['D:\NeuropixelData\CellReportsDataCode\Data\NeuropixelData\NeuronFiles\', session, '_1_sparse.mat'])
    bin_width = 0.05; % 1ms
    step = 0.05; % 1ms
    t_range = [-1 4];
    statecode_threshold = 7;
    bin_centers = t_range(1)+(step/2):step:t_range(2);

    [psth_mua, ~, target_trials] = compute_psth_sparse(MatData, ...
            bin_width, step, t_range, statecode_threshold, 'apply_stationary', false);
    %
    MatData_mua = MatData;
    for i = 1:numel(MatData_mua.trials)
        MatData_mua.trials(i).ss = MatData_mua.trials(i).ss(:, ~MatData_mua.single_units);
    end
    MatData_mua.stationary = MatData_mua.stationary(~MatData_mua.single_units, :);
    [psth_mua_plot, ~, ~] = compute_psth_sparse(MatData_mua, ...
        bin_width, step, t_range, statecode_threshold, 'apply_stationary', false);
    [psth_onlymua, ~, ~] = compute_psth_sparse(MatData_mua, ...
        bin_width/50, step/50, t_range, statecode_threshold, 'apply_stationary', false);

    MatData_single = MatData;
    for i = 1:numel(MatData_single.trials)
        MatData_single.trials(i).ss = MatData_single.trials(i).ss(:, MatData_single.single_units);
    end
    MatData_single.stationary = MatData_single.stationary(MatData_single.single_units, :);
    [psth_single_plot, ~, ~] = compute_psth_sparse(MatData_single, ...
        bin_width, step, t_range, statecode_threshold, 'apply_stationary', false);
    [psth_onlysingle, ~, ~] = compute_psth_sparse(MatData_single, ...
        bin_width/50, step/50, t_range, statecode_threshold, 'apply_stationary', false);

    if strcmpi(neuron_type, 'delay-unit')
        epoch = 'Delay';
        [~,~,SigNeurons] = xlsread(['D:\NeuropixelData\', session(1:6), '.xlsx'],[epoch, 'SigNeurons']);
        cids = str2double(SigNeurons(:,2));
        
        MatData_delay = MatData;
        c_ind = find(ismember(MatData_delay.cids, cids));
        
        for i = 1:numel(MatData_delay.trials)
            MatData_delay.trials(i).ss = MatData_delay.trials(i).ss(:, c_ind);
        end
        MatData_delay.stationary = MatData_delay.stationary(c_ind, :);
        [psth_delay, t_centers, target_trials] = compute_psth_sparse(MatData_delay, ...
            bin_width, step, t_range, statecode_threshold, 'apply_stationary', false);
        psth_delay = psth_delay(target_trials, :, :);
    end
    psth_mua = psth_mua(target_trials, :, :);
    psth_mua_plot = psth_mua_plot(target_trials, :, :);
    psth_single_plot = psth_single_plot(target_trials, :, :);
    psth_onlymua = psth_onlymua(target_trials, :, :);
    psth_onlysingle = psth_onlysingle(target_trials, :, :);

    t = (t_range(1):step:t_range(2)-bin_width) + bin_width/2;
    t2 = (t_range(1):step/50:t_range(2)-bin_width/50) + (bin_width/50)/2;


    
    ntr=0; fig=figure('Name','PIC138_class8_trial9');
    for tr=71%3%65%71%1:size(psth_mua, 1)
        ntr= ntr + 1;
        class = MatData.trials(target_trials(tr)).Class;
        tem_mua = squeeze(psth_mua(tr, :, :));
        tem_mua_plot = squeeze(psth_mua_plot(tr, :, :));
        tem_single_plot = squeeze(psth_single_plot(tr, :, :));

        tem_onlymua = squeeze(psth_onlymua(tr, :, :));
        tem_onlysingle = squeeze(psth_onlysingle(tr, :, :));

        tem_onlymua = tem_onlymua>0;
        tem_onlysingle = tem_onlysingle>0;
        
        total_height = 0.7;
        % Heights of subplots proportional to number of neurons
        height_su = total_height * (size(tem_onlysingle, 1) / (size(tem_onlysingle, 1) + size(tem_onlymua, 1)));
        height_mu = total_height * (size(tem_onlymua, 1) / (size(tem_onlysingle, 1) + size(tem_onlymua, 1)));
        
        ax1 = subplot('Position', [0.1, 0.99 - height_su, 0.8, height_su]);
        hold on;

        plotrasters(tem_onlysingle, t2, 'k', 'su')

        ax2 = subplot('Position', [0.1, 0.99 - height_su - height_mu - 0.01, 0.8, height_mu]);
        hold on;
        plotrasters(tem_onlymua,t2, 'k', 'mu')

        ax3 = subplot('Position', [0.1, 0.1, 0.8, 0.15]);
        hold on;

        mean_trial = mean(tem_mua);
        plot(t, mean_trial, 'color', 'k', 'LineWidth',2); hold on

        line(t_range, [mean(mean_trial(1:20)) mean(mean_trial(1:20))], 'LineStyle', '--', 'color', 'k'); hold on
        line([0 0], [min(mean(tem_mua))-1 max(mean(tem_mua))+1],'Color','k','LineStyle','-', 'LineWidth', 1); hold on
        line([0.5 0.5], [min(mean(tem_mua))-1 max(mean(tem_mua))+1],'Color','k','LineStyle','-', 'LineWidth', 1); hold on
        line([3.5 3.5], [min(mean(tem_mua))-1 max(mean(tem_mua))+1],'Color','k','LineStyle','-', 'LineWidth', 1); hold on
        xlim([-1 4])
        xlabel('Time from cue onset (s)', 'FontSize',15)
        ylabel('Firing rate (sp/s)', 'FontSize',15)
        yticks([floor(min(mean(tem_mua))):5:ceil(max(mean(tem_mua)))])
        p_value = perm_test_single_psth(mean_trial(1, 1:20), mean_trial(31:90))
%         significant_bins = sig_bins(mean_trial(1, 1:20), mean_trial(31:90))
        significant_bins = sig_bins(mean(tem_mua(:, 1:20),2), tem_mua(:, 31:90))


%         [h, p, ci, stat] = ttest2(mean_trial(1, 1:20), mean_trial(31:90))
        ax = gca; % Get current axes
        ax.FontSize = 15; % Set font size for y-axis
        axis tight
%         title('PSTH')
        box off
        if mod(ntr, 8)==0
%             pause
%             clf
            figure
        end
    end
    set(fig, 'Position', [726, -4, 670, 1000]);
%     sgtitle('PIC138_class8_trial159', 'FontSize',15)
end


function plotrasters(psth,t, colour,type)
[neuron_indices, time_bins] = find(psth);
hold on;
for i = 1:length(neuron_indices)
    n = neuron_indices(i);
    tb = time_bins(i);
    line([tb tb], [n-0.99 n-0.01], 'Color', colour, 'LineWidth', 0.5);
end

line([find(t>0,1) find(t>0,1)],[0 size(psth, 1)],'Color','k','LineStyle','-', 'LineWidth', 1)
line([find(t>0.5,1) find(t>0.5,1)],[0 size(psth, 1)],'Color','k','LineStyle','-', 'LineWidth', 1)
line([find(t>3.5,1) find(t>3.5,1)],[0 size(psth, 1)],'Color','k','LineStyle','-', 'LineWidth', 1)
xlim([-1 4])
if strcmpi(type, 'su')
    ylabel('Single-units', 'FontSize',15);
elseif strcmpi(type, 'mu')
    ylabel('Multi-units', 'FontSize',15);
end
yticks([10:30:size(psth, 1)])
ax = gca; % Get current axes
ax.XColor = 'none'; % Hide the x-axis

% Ensure the y-axis is visible and set its properties
ax.YColor = 'k'; % Set y-axis color to black (or any other desired color)
ax.FontSize = 15; % Set font size for y-axis
axis tight
hold on
end












