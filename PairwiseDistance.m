
 % find marker/cell type indices with threshold above eg 0.1 
 indices = find(abs(MPO)>0.1);MPO(indices)=[];
 
 % maintain only rows with these indices
 % insert XY coordinate - always the same
 MPO_XY = XY((indices),:); 
 
 % calculate all distance between cell pairs, eg
 D_Sox9_aSMA=pdist2(Sox9_XY,aSMA_XY); 
 D_Sox9_CALR=pdist2(Sox9_XY,CALR_XY);
 D_Sox9_CC3=pdist2(Sox9_XY,CC3_XY);
 D_Sox9_CD11c=pdist2(Sox9_XY,CD11c_XY);
 D_Sox9_Ly6G=pdist2(Sox9_XY,Ly6G_XY);
 D_Sox9_Arg1=pdist2(Sox9_XY,Arg1_XY);
 
 % plot in histogram form 
 figure (1);

%subplot (3,2,1);
histogram(D_CC3_Ly6G);
hold on 
xlim ([0 2000]);
ylim ([0 3200]);
title ('CC3-Ly6G', 'Fontsize', 12);