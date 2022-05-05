% data to present has a matrix fomat. collumns are biomarkers while rows
% are replicates split to two. If n=3 then first three rows (1:3) are controls
% and second 3 rows (4:6) are experimental condition

% average
ave (1, :) = mean (data(1:3, :)); % first ave row is the average of first three rows, n=3 for control
ave (2, :) = mean (data(4:6, :)); % second ave row is the average of the 4 through 6 row, n=3 for experimental. each column is a biomarker


%error_bar = std(data); % from the original code
error_bar (1, :) = std (data(1:3, :))/sqrt(3); % s.e. for n=3
error_bar (2, :) = std (data(4:6, :))/sqrt(3);

figure (2)
hold on
hb = bar(1:44,ave', 1, 'LineWidth',1.5); % number of markers and cell types presented, adapt for 2a (9), 3a, 4a (12)

hold on;
scatter(scatterDots(:,1), scatterDots(:,2), 'd', 'MarkerEdgeColor', [0.9290 0.6940 0.1250], 'linewidth', 2); %'MarkerEdgeColor', [0.55 0.55 0.55]

% For each set of bars, find the centers of the bars, and write error bars
pause(0.1); %pause allows the figure to be created
for ib = 1:numel(hb)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = hb(ib).XData+hb(ib).XOffset;
    errorbar(xData,ave(ib,:),error_bar(ib,:), 'k.','LineWidth',1.5)
end
ylim ([0 0.7])
%title('chemotaxis of CTRL MSCs', 'fontsize',14);
ylabel ('Rate of positive cells', 'fontsize',20)
set(gca, 'XTickLabel', {'-','CD11b+','-','CD45+','-', 'F4/80+','-','Ly6G+', '-', 'MHCII+', }, 'XTickLabelRotation', 35,'LineWidth',1.5, 'fontsize',25); %'XTickLabelRotation', 35,
legend ('random intratumoral','panobinostat ROI', 'Location','northeast','fontsize',30) 

%hold on;
%scatter(scatterDots(:,1), scatterDots(:,2), 'd', 'MarkerEdgeColor', [0.9290 0.6940 0.1250], 'linewidth', 2);
