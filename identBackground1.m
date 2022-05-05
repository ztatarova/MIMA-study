% code to identify if there are cells positive for all markers

% load or insert source data meaning 1) "fcsData" (cell ID, shape perimeter,
% markers, XY coordinates & 2) "percentage" of positive cells in the area

% extract marker (M) rage from fcsData only
% fcsData is numberX39 double
% percentage is 35x1 single, rate in fact
% no zero allowed in percentage file!

M = fcsData(:, 3:37);

% identify thershold based on the percentage of positivity

n=35; % number of markers
fcsThreshold = zeros (1,n); % allocate empty space for vecotr/matrix to speed up the for loop code

for i=1:n
    fcsThreshold(i)= min(maxk(M(:,i),ceil(size(M(:,i),1)*percentage(i,1))));
end

%% place 0 to all units in which expression is below the threshold value

idx = find(M < fcsThreshold);
M(sub2ind(size(M), idx)) = 0;

%% insert the "markers" matrix M back to the fcsData (create fcsDataB) to have cell ID and XY
% coordinates in teh same file

fcsDataB = fcsData;
fcsDataB(:, 3:37) = M;


%% finally identify and delete rows (each row is a cell) that are positive for all
% markers. This would very likely mean that we are looking at background
% signal. There is no cell that is positive. Again do this through logical
% indexing

id = all(fcsDataB > 0, 2);
fcsDataB (id, :) = [];

% this fcsDataB file is now ready for downstream data processing

