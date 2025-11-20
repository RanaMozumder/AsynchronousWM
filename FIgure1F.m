%% for individual session
clear all; close all; clc;
num_colors = 8;
hues = (0:num_colors-1) * (360 / num_colors); % 45 degrees apart
saturation = 0.75; % Reduced saturation for better aesthetics
value = 0.9; % Slightly reduced brightness
colors_rgb = zeros(num_colors, 3);
for i = 1:num_colors
    hsv = [hues(i) / 360, saturation, value];
    colors_rgb(i, :) = hsv2rgb(hsv);
end
epoch = 'Delay';
sessions = readtable("D:\NeuropixelData\CellReportsDataCode\Data\NeuropixelData\SessionInfo.xlsx", Sheet='EligibleSessions'); % reading sessions from the directory
for ses = 1%:height(sessions)
    raw = readtable(['D:\NeuropixelData\CellReportsDataCode\Data\NeuropixelData\SelectiveNeurons\', cell2mat(sessions{ses,1}), '.xlsx'],Sheet='cuedelay_selective');
    neuron = 0;
    for i = 2:height(raw)
        neuron = neuron + 1;
        filename = [cell2mat(raw.NeuronFileName(i)), '.mat'];
        Cue_pref(neuron)= raw.cue_pref(i);
        DIR='D:\NeuropixelData\single_neuron_files\PIC_VIK\';
        [~, ~, ~, SpkCount(neuron,:)] = ...
            SpikeCountVsClass([DIR, filename], epoch);
        [~, T(neuron,:), theta_pref(neuron,:)] = CalculateT(SpkCount(neuron,:), [0 45 90 135 180 225 270 315]);
    end
%     subplot(3, 4, ses-1)
    h= compass(T.*cos(deg2rad(theta_pref)), T.*sin(deg2rad(theta_pref)));
    for  n=1:size(T, 1)
        set(h(n), 'Color',colors_rgb(Cue_pref(n), :),'LineWidth', 1.5, 'MarkerSize', 2);
    end
    title(cell2mat(sessions{ses,1}))
    
    clear T_vec T theta_pref AlignedDelaySEM AlignedClasses AlignedDelayRates SpkCount raw Cue_pref
end


