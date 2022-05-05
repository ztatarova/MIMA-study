% marker presentation in XY coordinate system

% after defining threshold - now in fcsThreshold (1x29)
% plot marker expresion distribution and position in the selected area

% extract XY coordinates from fcsData file
XY=fcsData(:, 38:39);

%% Extract a marker from fcsData file
Arg1=fcsData(:, 3);

% old code
% identify threshold value by % of positive cells from fcsThreshold file
% Arg1Positive=maxk(Arg1,ceil(size(Arg1,1)*fcsThreshold(1,:))); % old
% fcsThreshold(3,2)=min(maxk(Arg1,ceil(size(Arg1,1)*fcsThreshold(3,1))));
% old code end

subplot (1,2,1);
hist(Arg1,200);
hold on;
line([fcsThreshold(1), fcsThreshold(1)], ylim, 'LineWidth', 2, 'Color', 'r');
hold on 
title ("Arg1 expresion and threshold for positivity", 'Fontsize', 14);
hold on

% identify indices with values above threshold
indices = find(abs(Arg1)>fcsThreshold(1));Arg1(indices)=[];
Arg1_XY = XY((indices),:); 
% Plot frequency of expresion rate and scatter plot to show if the cells are present on a correct spot

subplot (1,2,2);
scatter(-Arg1_XY(:,1), -Arg1_XY(:,2));
hold on 
title ("Arg1 position in ROI", 'Fontsize', 14);
hold on 
%xlim ([-5000 0]);
%ylim ([-4000 0]);

%% aSMA
aSMA=fcsData(:, 4);

subplot (1,2,1);
hist(aSMA,200);
hold on;
line([fcsThreshold(2), fcsThreshold(2)], ylim, 'LineWidth', 2, 'Color', 'r');
hold on 
title ("aSMA expresion and threshold for positivity", 'Fontsize', 14);
hold on

% identify indices with values above threshold
indices = find(abs(aSMA)>fcsThreshold(2));aSMA(indices)=[];
aSMA_XY = XY((indices),:); 
% Plot frequency of expresion rate and scatter plot to show if the cells are present on a correct spot

subplot (1,2,2);
scatter(-aSMA_XY(:,1), -aSMA_XY(:,2), 'filled');
hold on 
title ("aSMA position in ROI", 'Fontsize', 14);
hold on 
xlim ([-5000 0]);
ylim ([-5000 0]);

%% apply for any marker of interest

%% Sox9
Sox9=fcsData(:, 37);

subplot (1,2,1);
hist(Sox9,200);
hold on;
line([fcsThreshold(35), fcsThreshold(35)], ylim, 'LineWidth', 2, 'Color', 'r');
hold on 
title ("Sox9 expresion and threshold for positivity", 'Fontsize', 14);
hold on

% identify indices with values above threshold
indices = find(abs(Sox9)>fcsThreshold(35));Sox9(indices)=[];
Sox9_XY = XY((indices),:); 
% Plot frequency of expresion rate and scatter plot to show if the cells are present on a correct spot

subplot (1,2,2);
scatter(-Sox9_XY(:,1), -Sox9_XY(:,2));
hold on 
title ("Sox9 position in ROI", 'Fontsize', 14);

hold on 
xlim ([-5000 0]);
ylim ([-4000 0]);