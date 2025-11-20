function pev = pev_calculation(MatData)
% clc; close all;
% clearvars -except MatData

% df, degree of freedom = number of classes - 1
df = length(MatData.class) - 1;

% for cancelling biases, we consider even number of trials available for 
% each class
max_trials = min(cellfun(@(x) [length(x)], {MatData.class.ntr}));

% using a window size of 100 ms with 50 ms overlap
win_size = 0.1;
step_size = 0.05;

%bins
edge1 = -1: step_size: 5-win_size;
edge2 = -1+ win_size: step_size: 5;

% find x bar for all trials
total_trial = 0; 
psth_trial = nan(max_trials*(df+1), length(edge1));
class_trial = nan(max_trials*(df+1), length(edge1));
if ~isempty(MatData.class)
    % calculate
    for cl = 1:length(MatData.class)
        trials = randperm(length(MatData.class(cl).ntr));
        ntr = 0;
        for tr = trials(1:max_trials) %length(MatData.class(cl).ntr)
            ntr = ntr + 1; total_trial = total_trial + 1;
            TS = MatData.class(cl).ntr(tr).TS - MatData.class(cl).ntr(tr).Cue_onT;
            for bin=1:length(edge1)
                psth(cl, ntr, bin) = length(find(TS>edge1(bin) & TS<edge2(bin)))./win_size;
                psth_trial(total_trial, bin) = psth(cl, ntr, bin);
                class_trial(total_trial, bin) = cl;
            end
            clear TS
        end
    end
else
    disp('Empty file!')
end

% x_bar = squeeze(mean(psth, [1, 2]))';

%%
for bin = 1:length(edge1)
    % Perform ANOVA and get ANOVA table
    [p, tbl, anovaTable] = anova1(psth_trial(:, bin), class_trial(:, bin), "off");
    
    % Extract SS Between Groups, df, SS Total, and MSE
    SSBetweenGroups = tbl{2, 2};
    df = tbl{2, 3};
    SSTotal = tbl{4, 2};
    MSE = tbl{3, 4};
    
    % Calculate omega squared
    omegaSquared = (SSBetweenGroups - df * MSE) / (SSTotal + MSE);
%     SS_total = 0;
%     for j = 1:total_trial
%         SS_total = SS_total + (psth_trial(j,bin) - x_bar(1, bin)).^2;
%     end
%     SS_bg = 0;
%     for i=1:size(psth, 1)
%         SS_bg = SS_bg + max_trials.*((mean(psth(i, :,bin)) - x_bar(1, bin)).^2);
%     end
%     MSE = 0;
%     for i=1:size(psth, 1)
%         for j = 1:size(psth, 2)
%             MSE = MSE + (squeeze(psth(i, j, bin)) - mean(psth(i, :,bin))).^2;
%         end
%     end
    
%     pev(1, bin) = ((SS_bg - (df*MSE))./(SS_total + MSE));
    pev(1, bin) = omegaSquared;
    clear omegaSquared MSE SSTotal SSBetweenGroups df anovaTable
end

end












