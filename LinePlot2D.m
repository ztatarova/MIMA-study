% survival graphs

plot (x,y, 'Linewidth', 4);
set (gca, 'Linewidth', 1.5, 'fontsize',18, 'linewidth',1.5, 'box', 'off');
xlim ([0 70]);
ylim ([0 105]);

hold on
%legend({'control','venetoclax', 'venetoclax aPD1', 'venetoclax aCD40'}, 'Location', 'southwest', 'Fontsize', 14);
%title ('Immunization study, EMT6 model', 'Fontsize', 27);
xlabel ('Days since treatment initiation', 'Fontsize', 22);
ylabel ('% of surviving mice', 'Fontsize', 22)
