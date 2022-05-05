% calculate mean, standard deviation and standard error of data, collumns are conditions
ave = mean (data); 
sd = std (data);
se = sd/sqrt(4); % by n


figure
bar (ave, 0.75); % plot average
hold on;
sz = 70; % define size of each circle in the scatter
scatter(scatterDots(:,1), scatterDots(:,2), sz, 'd', 'MarkerEdgeColor', [0.9290 0.6940 0.1250], 'linewidth', 2);
hold on;
errorbar (ave, se,'k.','linewidth', 1.5);  % add standard error; 'k.' not connecting the errorbars
set(gca, 'XTickLabel', {'PEG control', 'panobinostat', 'panobinostat aGal3', 'panobinostat aLy6G', 'aGal3', 'aLy6G'}, 'fontsize', 20,'linewidth',1.5, 'XTickLabelRotation', 35, 'box', 'off');
ylabel('Mean intensity (px value)', 'fontsize',24);
%ylim ([0 0.55])
hold on;
% 'd' for diamond markers
%scatter(scatterDots(:,1), scatterDots(:,2), 'd', 'linewidth', 2); %
% alternative
%sz = 80; % define size of each circle in the scatter
%scatter(scatterDots(:,1), scatterDots(:,2), sz, 'k','linewidth', 2);