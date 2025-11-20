clc; clear all; close all;

% import session names from --> NeuropixelData\SessionInfo.xlsx
sessions = readtable('D:\NeuropixelData\CellReportsDataCode\Data\NeuropixelData\SessionInfo.xlsx', 'Sheet', 'EligibleSessions');
p=0; v=0;
for ses=1:height(sessions)
    % load session's neurons
    ses_name = cell2mat(sessions{ses, 1});
    load(['D:\NeuropixelData\beh\', ses_name, '_1.mat'])
    correct_trials = sum([AllData.trials.Statecode] == 7);
    total_trials = sum([AllData.trials.Statecode] > 4);
    if strcmpi(ses_name(1:3), 'PIC')
        p = p+1;
        pic(p) = correct_trials * 100 / total_trials;
    else
        v=v+1;
        vik(v) = correct_trials * 100 / total_trials;
    end
end
%% plot performance
figure('Position', [680, 558, 400, 600]); hold on
scatter(ones(length(pic)), pic, 'Marker', 'o', 'MarkerEdgeColor','k', 'SizeData', 150, 'LineWidth',2)
scatter(1.5*ones(length(vik)), vik, 'Marker', 'o', 'MarkerEdgeColor','k', 'SizeData', 150, 'LineWidth',2)
plot(0.8:0.1:1.2,mean(pic)*ones(5), 'k-', 'LineWidth', 2)
plot(1.3:0.1:1.7,mean(vik)*ones(5), 'k-', 'LineWidth', 2)
% boxplot([pic; vik], [ones(length(pic), 1); 2*ones(length(vik), 1)], 'Colors', 'k')
ylabel('Percentage correct')
xlim([0.7 1.8])
ylim([50 100])
xticks([1 1.5])
yticks([50 60 70 80 90 100])
xticklabels({'monkey P', 'monkey V'});
% set(x, 'FontSize', 16)
box off
title('Behavior Performance','FontSize',20)
set(gca, 'FontSize', 16)



