close all; clear all; clc;
epoch = 'cuedelay';
sessions = readtable('D:\NeuropixelData\CellReportsDataCode\Data\NeuropixelData\SessionInfo.xlsx', Sheet='EligibleSessions'); % reading sessions from the directory
win_size = 0.1;
step_size = 0.05;

%bins
edge1 = -1: step_size: 5-win_size;
edge2 = -1+ win_size: step_size: 5;
time = edge1+ win_size/2;

plots = []; ses_name = []; n=0;

figure
dark_colors = {
    '#8B0000', '#006400', '#00008B', '#8B008B', ...
    '#008B8B', '#FF8C00', '#9400D3', '#654321', ...
    '#2F4F4F', '#556B2F', '#483D8B', '#00CED1', ...
    '#2E2E2E'};


% Convert hex codes to RGB
dark_rgb = cellfun(@(x) hex2rgb(x), dark_colors, 'UniformOutput', false);
light_rgb = cellfun(@(c) blend_with_white(c, 0.7), dark_rgb, 'UniformOutput', false);
all_pev = {};
for ses = 1:height(sessions)
    if sessions{ses, 2}>80
        raw = readtable(['D:\NeuropixelData\CellReportsDataCode\Data\NeuropixelData\SelectiveNeurons\', cell2mat(sessions{ses,1}), '.xlsx'],Sheet=[epoch, '_selective']);
        if ~strcmpi(sessions{ses,1}, 'VIK093')
            n =n+1;
            session = cell2mat(sessions{ses,1})
            ses_name {n} = session;
            neuron = 0; pev = [];
            for i = 1:height(raw)
%                 i
                neuron = neuron + 1;
                filename = [cell2mat(raw.NeuronFileName(i)), '.mat'];
                DIR='D:\NeuropixelData\single_neuron_files\PIC_VIK\';
                load([DIR, filename])
                pev(neuron, :) = pev_calculation(MatData);
            end

            all_pev(ses) = {pev};
            if ~isempty(pev)
                pev_mean = mean(pev, 'omitnan')*100;
                pev_sem = (std(pev, 'omitnan')./sqrt(mean(sum(~isnan(pev)))))*100;
                shadedplot2(time, pev_mean+pev_sem, pev_mean-pev_sem, light_rgb{n}, light_rgb{n}); hold on
                p = plot(time, pev_mean, LineWidth=2,Color=dark_rgb{n}); hold on;
                plots = [plots p];

            end
        end
        clearvars -except ses sessions epoch time ses_name plots n dark_rgb light_rgb all_pev
    end
end

line([-2 5], [0 0], 'color', 'k', 'linestyle', '--')
line([0 0], [-2 25], 'color', 'k', 'linestyle', '--')
line([0.5 0.5 ], [-2 25], 'color', 'k', 'linestyle', '--')
line([3.5 3.5 ], [-2 25], 'color', 'k', 'linestyle', '--')
ylim([-2 25])
xlim([-1 5])
xlabel('Time from cue onset (s)', 'FontSize',15)
ylabel('\omega^2 PEV (%)', 'FontSize',15)
box off
legend(plots, ses_name)



function rgb = hex2rgb(hex)
hex = char(hex);
rgb = [hex2dec(hex(2:3)), hex2dec(hex(4:5)), hex2dec(hex(6:7))] / 255;
end

function shadedplot2(x, y1, y2, fillColor, edgeColor)
fill([x, fliplr(x)], [y1, fliplr(y2)], fillColor, 'FaceAlpha', 0.3, 'EdgeColor', edgeColor, 'EdgeAlpha', 0.3);
end

function blended_rgb = blend_with_white(rgb, blend_factor)
white = [1, 1, 1];
blended_rgb = (1 - blend_factor) * rgb + blend_factor * white;
end

