clc; clear all; close all;

monkey = 'PIC';
session = 138;
session = sprintf('%03d', session);
l_names = {'single-units'};
pdf_filename = sprintf('D:\\NeuropixelData\\Figures\\%s%s.pdf', monkey, session);
% find data from --> Data\NeuropixelData\NeuronFiles
load(['D:\NeuropixelData\CellReportsDataCode\Data\NeuropixelData\NeuronFiles\', monkey, session,'_1_sparse.mat'])

if strcmpi(monkey, 'PIC')
    chmap = load('D:\OptroUnity\KS2.5-main\NP1010_001180.mat'); %PIC
elseif strcmpi(monkey, 'VIK')
    chmap = load('D:\OptroUnity\KS2.5-main\NP1015_001180.mat'); % VIK
end
yc = chmap.ycoords(1:384);
xc = chmap.xcoords(1:384);

figure('Units','inches', 'Position', [0, 0, 7, 10]);
hold on

% s = scatter(xc, yc, 'ko', 'filled');
% s.MarkerFaceAlpha = 0.3;
% s.MarkerEdgeAlpha = 0.3;



n_wave_pt = size(MatData.cluster_waveforms, 2);
wave_plot_ratio = [0.4, 0.7];
wave_plot_edge  = round(n_wave_pt * wave_plot_ratio);
wave_plot_idx   = wave_plot_edge(1):wave_plot_edge(2);
n_wave_plot = numel(wave_plot_idx);
waveform_pos0 = ([1:n_wave_plot] - (n_wave_plot/2))/n_wave_plot*12;
waveform_scaler = 200;

cids = find(MatData.single_units);
for j = cids
    lw = 1.2;
    %         if mat_cell{i}.sync_excit(j)
    %             color_ = colors(1, :);
    %         else
    %             color_ = colors(5, :);
    %             lw = 0.9;
    %         end
    plot(MatData.clusterXs(j) + waveform_pos0, MatData.clusterYs(j) + MatData.cluster_waveforms(j, wave_plot_idx) * waveform_scaler, 'Color', rand(1, 3), 'LineWidth', lw);
end
%     ll = [];
%     for i_l = 1:size(colors, 1)
%         ll(i_l) = plot(0, 0, '-', 'Color', colors(i_l, :), 'DisplayName', l_names{i_l}, 'LineWidth', 1.2);
%     end
% legend(ll, 'Location', 'bestoutside');
if length(unique(xc))==4
    channelrow_x = sort(xc(1:4));
elseif length(unique(xc))==2
    channelrow_x = sort(xc(1:2));
end

equationText = sprintf('N = %d units', length(cids));
x_center = mean(xlim);
y_bottom = min(ylim);

text(x_center, y_bottom, equationText, 'FontSize', 14, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

xticks(channelrow_x)
ylabel('Vertical distance from the deepest channel (\mum)', 'Interpreter', 'tex')
xlabel('Horizontal position on the probe (\mum)', 'Interpreter', 'tex')
title([monkey, session]);
% ylim([-100 2500])
set(gca, 'FontSize', 16)
% exportgraphics(gcf, pdf_filename,"Append",true, 'ContentType', 'image', 'Resolution', 300);
% close(gcf);