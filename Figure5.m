%% plot simultaneous and shuffled decoding together
close all; clear all; clc
sessions = readtable('D:\NeuropixelData\CellReportsDataCode\Data\NeuropixelData\SessionInfo.xlsx', Sheet='EligibleSessions');

for s=1:2%height(sessions)
    ses = cell2mat(sessions.Session(s))
    realDecoding = importdata(['D:\NeuropixelData\CellReportsDataCode\Data\NeuropixelData\SMvsPseudo\', ses, '_1_logistic_sm_cv.mat']);
    
    load(['D:\NeuropixelData\CellReportsDataCode\Data\NeuropixelData\SMvsPseudo\', ses, '_1_logistic_shuff.mat']);
    for shuff = 1: length(accuracy)
        pseudoDecoding(shuff, :, :) = accuracy{shuff};
    end
    
    realAvg = mean(realDecoding, 1);   % 1x120
    pseudoAvg = squeeze(mean(pseudoDecoding, 2)); % 100x120
    
    pseudoMean = mean(pseudoAvg, 1);   % 1x120
    pseudoStd = std(pseudoAvg, [], 1); % 1x120
    
    disp('Signrank test')
    [p, h, stats] = signrank(realAvg, pseudoMean);
    ses_p(s) = p;
    % Step 3: Compute z-scores
    zScores = (realAvg - pseudoMean) ./ pseudoStd;
    
    % Step 4 (Optional): Calculate p-values
    pValues = 2 * (1 - normcdf(abs(zScores))); % Two-tailed
    
    % Optional multiple comparison correction
    % Example: Bonferroni correction (simple but conservative)
    adjustedPValues = pValues * 120; 
    adjustedPValues(adjustedPValues > 1) = 1; % Cap at 1
    
    t = -0.9750:0.05:5;
    % Plotting (example)
    fig = figure;
    % subplot(2,1,1);
    fill([t, flip(t)], [pseudoMean + pseudoStd, fliplr(pseudoMean - pseudoStd)], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName','Pseudo +- std');
    hold on
    plot(t, realAvg, 'b', 'DisplayName', 'Real');
    hold on;
    plot(t, pseudoMean, 'r--', 'DisplayName', 'Pseudo Mean');
    significanceThreshold = 0.05;
    significantTimePoints = find(adjustedPValues < significanceThreshold);
    plot(t(significantTimePoints), 0.95*ones(size(significantTimePoints)), 'k.', 'MarkerSize', 16, ...
        'DisplayName', 'Significant (p < 0.05)');
    legend;
    xlabel('Time Points');
    ylabel('Decoding Accuracy');
    title('Decoding Comparison');

    xline(0, 'color', 'k', 'Linestyle', '--')
    xline(0.5, 'color', 'k', 'Linestyle', '--')
    xline(3.5, 'color', 'k', 'Linestyle', '--')
    yline(0.125, 'color', 'k', 'Linestyle', '--')
    ylim([0 1])
    yticks([0.25 0.5 0.75 1])
    
    % subplot(2,1,2);
    % plot(1:120, adjustedPValues, 'k');
    % xlabel('Time Points');
    % ylabel('Adjusted p-value');
    % title('Statistical Significance Over Time');
    % yline(0.05, 'r--', 'p=0.05');
    box off
    savefig(['D:\NeuropixelData\Figures\Classification_onevsall\SMvsPS\Stat\', ses, '_smVSps_stat'])
    fclose all
    clearvars -except sessions s ses_p
end


%% quantify difference between simultaneous and shuffled decoding
% close all; 
clearvars; clc;
sessions = readtable('D:\NeuropixelData\CellReportsDataCode\Data\NeuropixelData\SessionInfo.xlsx', Sheet='EligibleSessions');

t = -0.9750:0.05:5;
tp_fix = find(t<0);
tp_cue = find(t>0 & t<0.5);
tp_delay = find(t>0.5 & t<3.5);
tp_saccade = find(t>3.5 & t<5);
for s=1:13
    ses = cell2mat(sessions.Session(s))
    % ses = 'PIC139_1';
    realDecoding = importdata(['D:\NeuropixelData\CellReportsDataCode\Data\NeuropixelData\SMvsPseudo\', ses, '_1_logistic_sm_cv.mat']);
    
    load(['D:\NeuropixelData\CellReportsDataCode\Data\NeuropixelData\SMvsPseudo\', ses, '_1_logistic_shuff.mat']);
    for shuff = 1: length(accuracy)
        pseudoDecoding(shuff, :, :) = accuracy{shuff};
    end
    
    
    % Assume realDecoding (5x120) and pseudoDecoding (100x5x120)
    
    % Step 1: Average across folds
    realAvg = mean(realDecoding, 1);   % 1x120
    pseudoAvg = squeeze(mean(pseudoDecoding, [2, 1]))'; % 100x120
    
    %difference at different epochs
    %fixation
%     diff_fix(s) = mean(sqrt((realAvg(1, tp_fix)-pseudoAvg(1,tp_fix)).^2));
    diff_fix(s) = mean(pseudoAvg(1, tp_fix)-realAvg(1,tp_fix));

    %cue
%     diff_cue(s) = mean(sqrt((realAvg(1, tp_cue)-pseudoAvg(1,tp_cue)).^2));
    diff_cue(s) = mean(pseudoAvg(1, tp_cue)-realAvg(1,tp_cue));

    %delay
%     diff_delay(s) = mean(sqrt((realAvg(1, tp_delay)-pseudoAvg(1,tp_delay)).^2));
    diff_delay(s) = mean(pseudoAvg(1, tp_delay)-realAvg(1,tp_delay));

    %saccade
%     diff_saccade(s) = mean(sqrt((realAvg(1, tp_saccade)-pseudoAvg(1,tp_saccade)).^2));
    diff_saccade(s) = mean(pseudoAvg(1, tp_saccade)-realAvg(1,tp_saccade));


end


%% box plots
% Combine data into a matrix
data = [diff_fix; diff_cue; diff_delay; diff_saccade]';  % Transpose to get 13x4 for boxchart

% Define group labels
groups = repelem(1:4, 13)';  % Each group has 13 values

% Create box plot
figure;
hold on;
boxchart(groups, data(:), 'BoxFaceColor', 'b', 'MarkerStyle', 'none');  % Box plot without default markers

% Overlay scatter plots
scatter(groups, data(:), 30, 'r', 'jitter', 'on', 'jitterAmount', 0.1);

% Customize plot
ylim([-0.07 0.06])
% yticks([0.02 0.04])
xticks(1:4);
xticklabels({'Fixation', 'Cue', 'Delay', 'Saccade'});
ylabel({'Difference in decoding accuracy', '(pseudo - simultaneous)'});
% title('Box Plots with Scatter Points');
set(gca, 'FontSize', 16)
%%
disp('fix vs cue')
[h, p, ci, stats] = ttest(diff_fix, diff_cue)
disp('fix vs delay')
[h, p, ci, stats] = ttest(diff_fix, diff_delay)
disp('fix vs saccade')
[h, p, ci, stats] = ttest(diff_fix, diff_saccade)

disp('cue vs delay')
[h, p, ci, stats] = ttest(diff_cue, diff_delay)
disp('cue vs saccade')
[h, p, ci, stats] = ttest(diff_cue, diff_saccade)

disp('saccade vs delay')
[h, p, ci, stats] = ttest(diff_delay, diff_saccade)
%%
plot([1,2], [-0.04, -0.04], 'k-')
scatter(1.5,-0.045, 'black','Marker','*')

plot([1,3], [-0.05, -0.05], 'k-')
scatter(2,-0.055, 'black','Marker','*')

plot([1,4], [-0.06, -0.06], 'k-')
scatter(2.5,-0.065, 'black','Marker','*')


