%profile plot with moving average filter to smooth the signal

% enter variables aSMA, Gal3 and k needed. k200 for 200px equals 100um

aSMAmean=mean(aSMA'); %average in cross section 
Gal3mean=mean(Gal3');
 
aSMAmovmean = movmean(aSMAmean, k); % moving average with 100um window
Gal3movmean = movmean(Gal3mean, k);

plot (Gal3movmean, 'm', 'LineWidth', 3); 
hold on; 
plot (aSMAmovmean,'g', 'LineWidth', 3 );
hold on;
legend ({'Gal3', 'aSMA'}, 'box','off');
set (gca, 'Linewidth', 1.5, 'fontsize',20, 'box','off');
xlim([850 3200]);