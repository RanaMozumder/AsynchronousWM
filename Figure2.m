%
clc; clearvars 
load('D:\NeuropixelData\CellReportsDataCode\Data\NeuropixelData\CrossTemporalDecoding\PIC138_1_logistic_wd.mat')
cmap = [...
    0 0 0;      % Black
    1 0 0;      % Red
    1 1 0];     % Yellow

% Interpolate the colormap to have smooth transitions
num_colors = 256;  % Number of colors in the colormap
custom_cmap = interp1(linspace(0, 1, size(cmap, 1)), cmap, linspace(0, 1, num_colors));


dur = [-1 5];
bin_width = 0.1; %100 ms
step = 0.05; %50 ms
t1 = dur(1):step:dur(2)-bin_width;
t = t1+bin_width/2;

figure
for tr=1:length(alpha_validate_cv)
    predictions2(tr, :, :) = double(squeeze(predictedLabels(tr,:,:))==alpha_validate_cv(tr));
end
avg2 = squeeze(mean(predictions2, 1));
subplot(3,3,5)
imagesc(t, t, flipud(avg2)); 
% colormap(cl_map)
% colormap(hot)
colormap(custom_cmap);

caxis([0.125 1])
% clim([0.5 1])

yticks([-0.5 2.5 3]);
yticklabels(flip([0 0.5 3.5]))
xticks([0 0.5 3.5]);
xticklabels([0 0.5 3.5])
xline(0, 'color', [1 1 1], 'Linestyle', '--')
xline(0.5, 'color', [1 1 1], 'Linestyle', '--')
xline(3.5, 'color', [1 1 1], 'Linestyle', '--')
yline(-0.5, 'color', [1 1 1], 'Linestyle', '--')
yline(3, 'color', [1 1 1], 'Linestyle', '--')
yline(2.5, 'color', [1 1 1], 'Linestyle', '--')
title('Average for across class')
xlabel('Testing time from cue (s)')
ylabel('Training time from cue (s)')
hold off
% clearvars -except MatData_single MatData custom_cmap

 
for class = 1:8
    ntr = 0;
    for tr=1:length(alpha_validate_cv)
        if alpha_validate_cv(1, tr) == class
            ntr = ntr + 1;
            predictions1(ntr, :, :) = double(squeeze(predictedLabels(tr,:,:))==alpha_validate_cv(tr));
        end
    end
    if ntr>0
    subplot(3,3,whichplot(class))
    avg1 = squeeze(mean(predictions1, 1));
    t = t;
    imagesc(t, t, flipud(avg1)); 
%     colormap(cl_map); 
%     colormap(hot)
    colormap(custom_cmap);

    caxis([0.125 1])
%     clim([0.5 1])

    yticks([-0.5 2.5 3]);
    yticklabels(flip([0 0.5 3.5]))
    xticks([0 0.5 3.5]);
    xticklabels([0 0.5 3.5])
    xline(0, 'color', [1 1 1], 'Linestyle', '--')
xline(0.5, 'color', [1 1 1], 'Linestyle', '--')
xline(3.5, 'color', [1 1 1], 'Linestyle', '--')
yline(-0.5, 'color', [1 1 1], 'Linestyle', '--')
yline(3, 'color', [1 1 1], 'Linestyle', '--')
yline(2.5, 'color', [1 1 1], 'Linestyle', '--')
%     title(sprintf('class=%d, #trials=%d', class, ntr))
    xlabel('Testing time from cue (s)')
    ylabel('Training time from cue (s)')
    end
%     clearvars -except MatData_single MatData class custom_cmap
    hold off
end

sgtitle('PIC138')

h = colorbar;

% Adjust the position of the colorbar to the right of the subplots
% Set the colorbar position relative to the figure
set(h, 'Position', [0.92, 0.1, 0.02, 0.8], 'ticks', [0.125 0.5 1]);



%%
function x = whichplot(class)
if class == 1
    x = 6;
elseif class == 2
    x = 3;
elseif class == 3
    x = 2;
elseif class == 4
    x = 1;
elseif class == 5
    x = 4;
elseif class == 6
    x = 7;
elseif class == 7
    x = 8;
elseif class == 8
    x = 9;
end
end









