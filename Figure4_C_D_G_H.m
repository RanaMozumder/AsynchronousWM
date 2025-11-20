% trained on PFC and tested on PFC and vice-versa
clc; clearvars
sessions = readtable('D:\NeuropixelData\CellReportsDataCode\Data\PFC_PPC_Data\BiasedODR.xlsx');
ss=[1:7];

ppc_on = 0; ppc_off = 0; pfc_on = 0; pfc_off = 0; off_lim = 3;
on_mat_pfc_all      = []; off_mat_pfc_all      = [];
on_mat_ppc_all      = []; off_mat_ppc_all      = [];
on_mat_pfc_null_all = []; off_mat_pfc_null_all = [];
on_mat_ppc_null_all = []; off_mat_ppc_null_all = [];

psth_neuron_pfc = []; psth_on_neuron_pfc = []; psth_off_neuron_pfc = []; 
psth_neuron_ppc = []; psth_on_neuron_ppc = []; psth_off_neuron_ppc = []; 
for s = ss
    ses = cell2mat(sessions.Ses(s))
    load(['D:\NeuropixelData\CellReportsDataCode\Data\PFC_PPC_Data\PFC_PPC_CTD\', ses(1:6),'_bc_logistic_shuffled.mat']);

    on_mat_pfc           = nan(size(predictedScores_pfc, 1), size(pfc_psth_train, 2), length(unique(alpha_validate_cv)));
    on_mat_null_pfc      = nan(1000, size(predictedScores_pfc, 1), size(pfc_psth_train, 2), length(unique(alpha_validate_cv)));
    off_mat_null_pfc     = nan(1000, size(predictedScores_pfc, 1), size(pfc_psth_train, 2), length(unique(alpha_validate_cv)));
    off_mat_pfc          = nan(size(predictedScores_pfc, 1), size(pfc_psth_train, 2), length(unique(alpha_validate_cv)));
    
    on_mat_ppc           = nan(size(predictedScores_ppc, 1), size(ppc_psth_train, 2), length(unique(alpha_validate_cv)));
    on_mat_null_ppc      = nan(1000, size(predictedScores_ppc, 1), size(ppc_psth_train, 2), length(unique(alpha_validate_cv)));
    off_mat_null_ppc     = nan(1000, size(predictedScores_ppc, 1), size(ppc_psth_train, 2), length(unique(alpha_validate_cv)));
    off_mat_ppc          = nan(size(predictedScores_ppc, 1), size(ppc_psth_train, 2), length(unique(alpha_validate_cv)));
    
    rng(0);
    for shuff=1:1000
    alpha_null(shuff,:) = alpha_validate_cv(randperm(length(alpha_validate_cv)));
    end
    for cl=1:length(unique(alpha))
        psth_held_ses_temp_pfc = []; psth_held_ses_temp_ppc = []; i=0;
        for tr=1:size(predictedScores_pfc, 1)
            if alpha_validate_cv(tr)==cl
                i=i+1;
                trial_n= tested_trial(tr);
                psth_held_ses_temp_pfc(i, :, :) = squeeze(psth_pfc(trial_n,:,:));
                psth_held_ses_temp_ppc(i, :, :) = squeeze(psth_ppc(trial_n,:,:));
            end
        end
        psth_held_on_pfc      = nan(size(psth_held_ses_temp_pfc));
        psth_held_off_pfc     = nan(size(psth_held_ses_temp_pfc));

        psth_held_on_ppc      = nan(size(psth_held_ses_temp_ppc));
        psth_held_off_ppc     = nan(size(psth_held_ses_temp_ppc));
        i=0;
        for tr=1:size(predictedScores_pfc, 1)
            if alpha_validate_cv(tr)==cl
                i=i+1;
                trial_n = tested_trial(tr);

                psth_heldout_pfc = squeeze(pfc_psth_train(trial_n,:,:));
                psth_heldout_ppc = squeeze(ppc_psth_train(trial_n,:,:));

                psth_heldout_pfc_og = squeeze(psth_pfc(trial_n,:,:));
                psth_heldout_ppc_og = squeeze(psth_ppc(trial_n,:,:));

                true_confidence_pfc = double(squeeze(predictedScores_pfc(tr, :, 2)));
                true_confidence_ppc = double(squeeze(predictedScores_ppc(tr, :, 2)));

                for iter=1:size(predictedScoresShuff_pfc, 3)
                    null_pfc(iter, :) = double(squeeze(predictedScoresShuff_pfc(tr, :, iter, 2)));
                    null_ppc(iter, :) = double(squeeze(predictedScoresShuff_ppc(tr, :, iter, 2)));
                end
                [pvals_on_pfc, ~, clustIdx_on_pfc, rawP_on_pfc, offz_pfc] = clusterMassOneSampZ([true_confidence_pfc;null_pfc]);
                [pvals_on_ppc, ~, clustIdx_on_ppc, rawP_on_ppc, offz_ppc] = clusterMassOneSampZ([true_confidence_ppc;null_ppc]);

                %%PFC
                %plot off states
                [clustIdx, nClusts] =  bwlabel(offz_pfc);
                sig_c = 0;  temp_c = [];
                for c = 1:nClusts
                    clusterTimePoints = t(clustIdx == c);
                    clusterTimePoints = clusterTimePoints(clusterTimePoints>=0.5 & clusterTimePoints<=3.5);
                    if numel(clusterTimePoints)>=off_lim
                        [~, indices]         = ismember(clusterTimePoints, t);
                        temp_c(:, sig_c + 1) = mean(psth_heldout_pfc(:, indices), 2);
                        psth_held_off_pfc(i, :, indices) = psth_heldout_pfc_og(:, indices);
                        sig_c                = sig_c +1;
                        pfc_off              = pfc_off + 1;
                        clear indices
                    end
                end
                if ~isempty(temp_c)
                    off_mat_pfc(tr, :, alpha_validate_cv(tr)) = mean(temp_c, 2);
                    for shuff =1:1000
                    off_mat_null_pfc(shuff, tr, :, alpha_null(shuff,tr))   = mean(temp_c, 2);
                    end
                end
                clear temp_c offz_pfc nClusts clustIdx sig_c clusterTimePoints

                %plot on states
                sigClusters = find(pvals_on_pfc < 0.05);
                sig_c = 0; temp_c = [];
                for c = 1:length(sigClusters)
                    clusterTimePoints = t(clustIdx_on_pfc == sigClusters(c)); % Time points in the cluster
                    clusterTimePoints = clusterTimePoints(clusterTimePoints>=0.5 & clusterTimePoints<=3.5);
                    if ~isempty(clusterTimePoints)
                        [~, indices]         = ismember(clusterTimePoints, t);
                        temp_c(:, sig_c + 1) = mean(psth_heldout_pfc(:, indices), 2);
                        psth_held_on_pfc(i, :, indices) = psth_heldout_pfc_og(:, indices);
                        sig_c                = sig_c +1;
                        pfc_on               = pfc_on + 1;
                        clear indices
                    end
                end
                if ~isempty(temp_c)
                    on_mat_pfc(tr, :, alpha_validate_cv(tr)) = mean(temp_c, 2);
                    for shuff =1:1000
                    on_mat_null_pfc(shuff, tr, :, alpha_null(shuff,tr))   = mean(temp_c, 2);
                    end
                end
                clear temp_c sigClusters sig_c clusterTimePoints

                %%PPC
                %plot off states
                [clustIdx, nClusts] =  bwlabel(offz_ppc);
                sig_c = 0;  temp_c = [];
                for c = 1:nClusts
                    clusterTimePoints = t(clustIdx == c);
                    clusterTimePoints = clusterTimePoints(clusterTimePoints>=0.5 & clusterTimePoints<=3.5);
                    if numel(clusterTimePoints)>=off_lim
                        [~, indices]         = ismember(clusterTimePoints, t);
                        temp_c(:, sig_c + 1) = mean(psth_heldout_ppc(:, indices), 2);
                        psth_held_off_ppc(i, :, indices) = psth_heldout_ppc_og(:, indices);
                        sig_c                = sig_c +1;
                        ppc_off              = ppc_off + 1;
                        clear indices
                    end
                end
                if ~isempty(temp_c)
                    off_mat_ppc(tr, :, alpha_validate_cv(tr)) = mean(temp_c, 2);
                    for shuff =1:1000
                    off_mat_null_ppc(shuff, tr, :, alpha_null(shuff,tr))   = mean(temp_c, 2);
                    end
                end
                clear temp_c offz_ppc nClusts clustIdx sig_c clusterTimePoints

                %plot on states
                sigClusters = find(pvals_on_ppc < 0.05);
                sig_c = 0;  temp_c = [];
                for c = 1:length(sigClusters)
                    clusterTimePoints = t(clustIdx_on_ppc == sigClusters(c));
                    clusterTimePoints = clusterTimePoints(clusterTimePoints>=0.5 & clusterTimePoints<=3.5);
                    if ~isempty(clusterTimePoints)
                        [~, indices]         = ismember(clusterTimePoints, t);
                        temp_c(:, sig_c + 1) = mean(psth_heldout_ppc(:, indices), 2);
                        psth_held_on_ppc(i, :, indices) = psth_heldout_ppc_og(:, indices);
                        sig_c                = sig_c +1;
                        ppc_on               = ppc_on + 1;
                        clear indices
                    end
                end
                if ~isempty(temp_c)
                    on_mat_ppc(tr, :, alpha_validate_cv(tr)) = mean(temp_c, 2);
                    for shuff =1:1000
                    on_mat_null_ppc(shuff, tr, :, alpha_null(shuff,tr))   = mean(temp_c, 2);
                    end
                end
                clear temp_c sigClusters sig_c clusterTimePoints

                psth_held_ses_temp_pfc(i, :, :) = psth_heldout_pfc_og;
                psth_held_ses_temp_ppc(i, :, :) = psth_heldout_ppc_og;
            end
        end
        psth_ses_pfc{cl}         = {squeeze(mean(psth_held_ses_temp_pfc, 1, 'omitnan'))};
        psth_on_ses_pfc{cl}      = {squeeze(mean(psth_held_on_pfc, 1, 'omitnan'))};
        psth_off_ses_pfc{cl}     = {squeeze(mean(psth_held_off_pfc, 1, 'omitnan'))};

        psth_ses_ppc{cl}         = {squeeze(mean(psth_held_ses_temp_ppc, 1, 'omitnan'))};
        psth_on_ses_ppc{cl}      = {squeeze(mean(psth_held_on_ppc, 1, 'omitnan'))};
        psth_off_ses_ppc{cl}     = {squeeze(mean(psth_held_off_ppc, 1, 'omitnan'))};
    end
    if size(on_mat_pfc, 2)>=80
        on_mat_pfc_all            = [on_mat_pfc_all;            squeeze(mean(on_mat_pfc,       1, 'omitnan'))];
        off_mat_pfc_all           = [off_mat_pfc_all;           squeeze(mean(off_mat_pfc,      1, 'omitnan'))];
        if isempty(on_mat_pfc_null_all)
            on_mat_pfc_null_all = squeeze(mean(on_mat_null_pfc,  2, 'omitnan'));
        else
            on_mat_pfc_null_all = cat(2, on_mat_pfc_null_all, squeeze(mean(on_mat_null_pfc,  2, 'omitnan')));
        end
        if isempty(off_mat_pfc_null_all)
            off_mat_pfc_null_all = squeeze(mean(off_mat_null_pfc,  2, 'omitnan'));
        else
            off_mat_pfc_null_all = cat(2, off_mat_pfc_null_all, squeeze(mean(off_mat_null_pfc,  2, 'omitnan')));
        end
        
        i=0;
        for n = 1: size(psth_ses_pfc{1, 1}{1, 1}, 1)
            i=i+1; temp = []; temp_on = []; temp_off = [];
            [~, maxind] = max(sum([on_mat_pfc_all(end-size(psth_ses_pfc{1, 1}{1, 1}, 1)+i,:);off_mat_pfc_all(end-size(psth_ses_pfc{1, 1}{1, 1}, 1)+i,:)], 1));
            
            temp         = psth_ses_pfc{1, maxind}{1, 1};
            temp_on      = psth_on_ses_pfc{1, maxind}{1, 1};
            temp_off     = psth_off_ses_pfc{1, maxind}{1, 1};

            psth_neuron_ses_pfc(n, :)        = temp(n, :);
            psth_on_neuron_ses_pfc(n, :)     = temp_on(n, :);
            psth_off_neuron_ses_pfc(n, :)    = temp_off(n, :);
        end
        psth_neuron_pfc         = [psth_neuron_pfc;         psth_neuron_ses_pfc];
        psth_on_neuron_pfc      = [psth_on_neuron_pfc;      psth_on_neuron_ses_pfc];
        psth_off_neuron_pfc     = [psth_off_neuron_pfc;     psth_off_neuron_ses_pfc];
    end
    if size(on_mat_ppc, 2)>=80
        on_mat_ppc_all            = [on_mat_ppc_all;            squeeze(mean(on_mat_ppc,       1, 'omitnan'))];
        off_mat_ppc_all           = [off_mat_ppc_all;           squeeze(mean(off_mat_ppc,      1, 'omitnan'))];
        if isempty(on_mat_ppc_null_all)
            on_mat_ppc_null_all = squeeze(mean(on_mat_null_ppc,  2, 'omitnan'));
        else
            on_mat_ppc_null_all = cat(2, on_mat_ppc_null_all, squeeze(mean(on_mat_null_ppc,  2, 'omitnan')));
        end
        if isempty(off_mat_ppc_null_all)
            off_mat_ppc_null_all = squeeze(mean(off_mat_null_ppc,  2, 'omitnan'));
        else
            off_mat_ppc_null_all = cat(2, off_mat_ppc_null_all, squeeze(mean(off_mat_null_ppc,  2, 'omitnan')));
        end
           
        i=0;
        for n = 1: size(psth_ses_ppc{1, 1}{1, 1}, 1)
            i=i+1; temp = []; temp_on = []; temp_off = []; 
            [~, maxind] = max(sum([on_mat_ppc_all(end-size(psth_ses_ppc{1, 1}{1, 1}, 1)+i,:);off_mat_ppc_all(end-size(psth_ses_ppc{1, 1}{1, 1}, 1)+i,:)], 1));
            
            temp         = psth_ses_ppc{1, maxind}{1, 1};
            temp_on      = psth_on_ses_ppc{1, maxind}{1, 1};
            temp_off     = psth_off_ses_ppc{1, maxind}{1, 1};

            psth_neuron_ses_ppc(n, :)     = temp(n, :);
            psth_on_neuron_ses_ppc(n, :)  = temp_on(n, :);
            psth_off_neuron_ses_ppc(n, :) = temp_off(n, :);
        end
        psth_neuron_ppc         = [psth_neuron_ppc;         psth_neuron_ses_ppc];
        psth_on_neuron_ppc      = [psth_on_neuron_ppc;      psth_on_neuron_ses_ppc];
        psth_off_neuron_ppc     = [psth_off_neuron_ppc;     psth_off_neuron_ses_ppc];
    end
    
    clearvars -except t ss off_mat_pfc_null_all on_mat_pfc_null_all off_mat_pfc_all on_mat_pfc_all ...
        off_mat_ppc_null_all on_mat_ppc_null_all off_mat_ppc_all on_mat_ppc_all ...
        ses pfc_on pfc_off ppc_on ppc_off s sessions off_lim psth_neuron_pfc psth_on_neuron_pfc psth_off_neuron_pfc ...
        psth_neuron_ppc psth_on_neuron_ppc psth_off_neuron_ppc psth_nothing_neuron_ppc psth_nothing_neuron_pfc ...
        nothing_mat_pfc_all nothing_mat_pfc_null_all avg_mat_pfc_all avg_mat_pfc_null_all ...
        nothing_mat_ppc_all nothing_mat_ppc_null_all avg_mat_ppc_all avg_mat_ppc_null_all
end
%%
t_start = find(t>0.5,1);
t_end   = find(t>3.5,1)-1;
figure; sgtitle('PFC')
subplot(2,1, 1)
[gauss_on_fit_pfc_pfc, gauss_off_fit_pfc_pfc, rotated_on_mat_pfc_pfc, rotated_off_mat_pfc_pfc] =...
    plottuning(on_mat_pfc_all, off_mat_pfc_all, on_mat_pfc_null_all, ...
    off_mat_pfc_null_all, pfc_on, pfc_off, [0.85,0.33,0.10], [0.00,0.45,0.74]);
title(sprintf('Tuning Function (n = %d)', size(on_mat_pfc_all, 1)))


subplot(2,1,2); hold on;
avg_psth = mean(psth_neuron_pfc, 1, 'omitnan');
std_psth = std(psth_neuron_pfc, 1, 'omitnan')./sqrt(size(psth_neuron_pfc, 1));

avg_psth_on = mean(psth_on_neuron_pfc, 1, 'omitnan');
std_psth_on = std(psth_on_neuron_pfc, 1, 'omitnan')./sqrt(size(psth_on_neuron_pfc, 1));

avg_psth_off = mean(psth_off_neuron_pfc, 1, 'omitnan');
std_psth_off = std(psth_off_neuron_pfc, 1, 'omitnan')./sqrt(size(psth_off_neuron_pfc, 1));

fill([t', flip(t')], [avg_psth + std_psth, fliplr(avg_psth - std_psth)], 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none');hold on
fill([t(t_start:t_end)', flip(t(t_start:t_end)')], [avg_psth_on(t_start:t_end) + std_psth_on(t_start:t_end), fliplr(avg_psth_on(t_start:t_end) - std_psth_on(t_start:t_end))], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');hold on
fill([t(t_start:t_end)', flip(t(t_start:t_end)')], [avg_psth_off(t_start:t_end) + std_psth_off(t_start:t_end), fliplr(avg_psth_off(t_start:t_end) - std_psth_off(t_start:t_end))], 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none');hold on

p1 = plot(t, avg_psth,         'Color','k'      , 'LineWidth', 1);
p2 = plot(t, avg_psth_on,      'Color', [0.85,0.33,0.10], 'LineWidth', 1);
p3 = plot(t, avg_psth_off,     'Color', [0.00,0.45,0.74], 'LineWidth', 1);

xlim([-0.5 3.5])
xline(0,   'linestyle', '--')
xline(0.5, 'linestyle', '--')
% xline(3.5, 'linestyle', '--')
% yline(0,   'linestyle', '--')
legend([p1, p2, p3], 'average', 'on', 'off')
title(sprintf('PSTH (n=%d, on=%d, off=%d)', size(on_mat_pfc_all, 1), pfc_on, pfc_off))

xlabel('Time from cue onset (s)')
ylabel('Firing rate (Hz)')

%PPC
figure; sgtitle('PPC')
subplot(2, 1, 1)
[gauss_on_fit_ppc_ppc, gauss_off_fit_ppc_ppc, rotated_on_mat_ppc_ppc, rotated_off_mat_ppc_ppc] = ...
    plottuning(on_mat_ppc_all, off_mat_ppc_all, ...
    on_mat_ppc_null_all, off_mat_ppc_null_all, ppc_on, ppc_off, [0.85,0.33,0.10], [0.00,0.45,0.74]);
title(sprintf('Tuning Function (n = %d)', size(on_mat_ppc_all, 1)))


subplot(2, 1, 2); hold on;
avg_psth = mean(psth_neuron_ppc, 1, 'omitnan');
std_psth = std(psth_neuron_ppc, 1, 'omitnan')./sqrt(size(psth_neuron_ppc, 1));

avg_psth_on = mean(psth_on_neuron_ppc, 1, 'omitnan');
std_psth_on = std(psth_on_neuron_ppc, 1, 'omitnan')./sqrt(size(psth_on_neuron_ppc, 1));

avg_psth_off = mean(psth_off_neuron_ppc, 1, 'omitnan');
std_psth_off = std(psth_off_neuron_ppc, 1, 'omitnan')./sqrt(size(psth_off_neuron_ppc, 1));

fill([t', flip(t')], [avg_psth + std_psth, fliplr(avg_psth - std_psth)], 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none');hold on
fill([t(t_start:t_end)', flip(t(t_start:t_end)')], [avg_psth_on(t_start:t_end) + std_psth_on(t_start:t_end), fliplr(avg_psth_on(t_start:t_end) - std_psth_on(t_start:t_end))], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');hold on
fill([t(t_start:t_end)', flip(t(t_start:t_end)')], [avg_psth_off(t_start:t_end) + std_psth_off(t_start:t_end), fliplr(avg_psth_off(t_start:t_end) - std_psth_off(t_start:t_end))], 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none');hold on

p1 = plot(t, avg_psth,         'Color','k'      , 'LineWidth', 1);
p2 = plot(t, avg_psth_on,      'Color', [0.85,0.33,0.10], 'LineWidth', 1);
p3 = plot(t, avg_psth_off,     'Color', [0.00,0.45,0.74], 'LineWidth', 1);

xlim([-0.5 3.5])
xline(0,   'linestyle', '--')
xline(0.5, 'linestyle', '--')
% xline(3.5, 'linestyle', '--')
% yline(0,   'linestyle', '--')
legend([p1, p2, p3], 'average', 'on', 'off')
% legend([p1, p2, p3, p4], 'average', 'on', 'off', 'neither')
title(sprintf('PSTH (n=%d, on=%d, off=%d)', size(on_mat_ppc_all, 1), ppc_on, ppc_off))

xlabel('Time from cue onset (s)')
ylabel('Firing rate (Hz)')

save('D:\NeuropixelData\Figures\OnOff\SameAreaFits_old.mat', 'gauss_on_fit_pfc_pfc', ...
    'gauss_off_fit_pfc_pfc','gauss_on_fit_ppc_ppc', 'gauss_off_fit_ppc_ppc')

%%
figure;
subplot 211
tun_depth_diff = (rotated_on_mat_pfc_pfc(:,5)-rotated_on_mat_pfc_pfc(:,1)) - ...
    (rotated_off_mat_pfc_pfc(:,5)-rotated_off_mat_pfc_pfc(:,1));
histogram(tun_depth_diff, 'BinWidth',0.1)
disp('PFC PFC tuning')
disp('signrank test')
[p, h, stats] = signrank(tun_depth_diff, 0)
disp('ttest test')
[h, p, ci, stats] = ttest(tun_depth_diff, 0)
title(sprintf('PFC PFC, p=%.3f', p))

subplot 212
tun_depth_diff = (rotated_on_mat_ppc_ppc(:,5)-rotated_on_mat_ppc_ppc(:,1)) - ...
    (rotated_off_mat_ppc_ppc(:,5)-rotated_off_mat_ppc_ppc(:,1));
histogram(tun_depth_diff, 'BinWidth',0.1)
disp('PPC PPC tuning')
disp('signrank test')
[p, h, stats] = signrank(tun_depth_diff, 0)
disp('ttest test')
[h, p, ci, stats] = ttest(tun_depth_diff, 0)
title(sprintf('PPC PPC, p=%.3f', p))
%%
clc
tun_diff = (rotated_on_mat_pfc_pfc(:,5)-rotated_off_mat_pfc_pfc(:,5));
disp('PFC PFC tuning')
disp('signrank test')
[p, h, stats] = signrank(tun_diff, 0)
disp('ttest test')
[h, p, ci, stats] = ttest(tun_diff, 0)


tun_diff = (rotated_on_mat_ppc_ppc(:,5)-rotated_off_mat_ppc_ppc(:,5));
disp('PPC PPC tuning')
disp('signrank test')
[p, h, stats] = signrank(tun_diff, 0)
disp('ttest test')
[h, p, ci, stats] = ttest(tun_diff, 0)

%% Stat testing
clc
%PFC
% On Vs Off
data1 = mean(psth_on_neuron_pfc(:, t_start:t_end), 2, 'omitnan');
data2 = mean(psth_off_neuron_pfc(:, t_start:t_end), 2, 'omitnan');

data1 = data1(~isnan(data1));
data2 = data2(~isnan(data2));

disp('On Vs Off')
[h, p, ci, stats] = ttest(data1, data2)    

% On Vs Avg

data1 = mean(psth_on_neuron_pfc(:, t_start:t_end), 2, 'omitnan');
data2 = mean(psth_neuron_pfc(:, t_start:t_end), 2, 'omitnan');

data1 = data1(~isnan(data1));
data2 = data2(~isnan(data2));

disp('On Vs Avg')
[h, p, ci, stats] = ttest(data1, data2)



% Avg Vs Off

data1 = mean(psth_neuron_pfc(:, t_start:t_end), 2, 'omitnan');
data2 = mean(psth_off_neuron_pfc(:, t_start:t_end), 2, 'omitnan');

data1 = data1(~isnan(data1));
data2 = data2(~isnan(data2));

disp('Avg Vs Off')
[h, p, ci, stats] = ttest(data1, data2)



%PPC
% On Vs Off

data1 = mean(psth_on_neuron_ppc(:, t_start:t_end), 2, 'omitnan');
data2 = mean(psth_off_neuron_ppc(:, t_start:t_end), 2, 'omitnan');

data1 = data1(~isnan(data1));
data2 = data2(~isnan(data2));

disp('On Vs Off')
[h, p, ci, stats] = ttest(data1, data2)    

% On Vs Avg

data1 = mean(psth_on_neuron_ppc(:, t_start:t_end), 2, 'omitnan');
data2 = mean(psth_neuron_ppc(:, t_start:t_end), 2, 'omitnan');

data1 = data1(~isnan(data1));
data2 = data2(~isnan(data2));

disp('On Vs Avg')
[h, p, ci, stats] = ttest(data1, data2)



% Avg Vs Off

data1 = mean(psth_neuron_ppc(:, t_start:t_end), 2, 'omitnan');
data2 = mean(psth_off_neuron_ppc(:, t_start:t_end), 2, 'omitnan');

data1 = data1(~isnan(data1));
data2 = data2(~isnan(data2));

disp('Avg Vs Off')
[h, p, ci, stats] = ttest(data1, data2)

%% Reduction
((max(gauss_on_fit_pfc_pfc)-min(gauss_on_fit_pfc_pfc))-...
    (max(gauss_off_fit_pfc_pfc)-min(gauss_off_fit_pfc_pfc)))/(max(gauss_on_fit_pfc_pfc)-min(gauss_on_fit_pfc_pfc))



((max(gauss_on_fit_ppc_ppc)-min(gauss_on_fit_ppc_ppc))-...
    (max(gauss_off_fit_ppc_ppc)-min(gauss_off_fit_ppc_ppc)))/(max(gauss_on_fit_ppc_ppc)-min(gauss_on_fit_ppc_ppc))


