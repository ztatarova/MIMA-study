% IntersectCellTypes3 is to find and display standard cell types (rather than
% markers) in XY coordinate system. 

% allows cell state classification. eg N1, N2 neutrophils.

% approach also potentiates identification of new biomarkers or response and resistance based on
% standard cell type classification. eg Galectin-3 as a biomarker enriching
% CSCs in BC; ICAM-1 preferably on immune or tumor cells

% Run IdentBackground1

% extract XY coordinates of all cells from fcsData file
% should be the last two collumns in the source csv array
XY=fcsData(:, 38:39);

%% Epcam negative

%extract collumn coresponding to Epcam
Epcam=fcsData(:, 18);

% old code
% identify threshold value by % of positive cells from fcsThreshold file
% Arg1Positive=maxk(Arg1,ceil(size(Arg1,1)*fcsThreshold(1,:))); % old
% fcsThreshold(3,2)=min(maxk(Arg1,ceil(size(Arg1,1)*fcsThreshold(3,1))));

% identify indices with values below threshold
indicesEpcamNeg = find(abs(Epcam)<fcsThreshold(16));Epcam(indicesEpcamNeg)=[];
EpcamNeg_XY = XY((indicesEpcamNeg),:); 

% plot to be sure
scatter(-EpcamNeg_XY(:,1), -EpcamNeg_XY(:,2),'filled');
hold on 
xlim ([-5000 0]);
ylim ([-4000 0]);

%% Epcam positive

% pre-delete Epcam variable for overlap option 
Epcam=fcsData(:, 18);

% identify indices with values below threshold
indicesEpcam = find(abs(Epcam)>fcsThreshold(16));Epcam(indicesEpcam)=[];
Epcam_XY = XY((indicesEpcam),:); 

% plot to be sure
scatter(-Epcam_XY(:,1), -Epcam_XY(:,2),'filled');
hold on 
xlim ([-5000 0]);
ylim ([-4000 0]);

%% CD45 negative

CD45=fcsData(:, 12);

% identify indices with values below threshold
indicesCD45Neg = find(abs(CD45)<fcsThreshold(10));CD45(indicesCD45Neg)=[];
CD45Neg_XY = XY((indicesCD45Neg),:); 

% Plot scatter plot to show if the cells are present on a correct spot
scatter(-CD45Neg_XY(:,1), -CD45Neg_XY(:,2),'filled');
hold on 
xlim ([-5000 0]);
ylim ([-4000 0]);

%% CD45 positive

% pre-delete CD45 variable
CD45=fcsData(:, 12);

% identify indices with values below threshold
indicesCD45 = find(abs(CD45)>fcsThreshold(10));CD45(indicesCD45)=[];
CD45_XY = XY((indicesCD45),:); 

% Plot scatter plot to show if the cells are present on a correct spot
scatter(-CD45_XY(:,1), -CD45_XY(:,2), 'filled', 'MarkerEdgeColor',[0.4667    0.6745    0.1882],'MarkerFaceColor',[0.4667    0.6745    0.1882] );
hold on 
xlim ([-5000 0]);
ylim ([-4000 0]);

%% CELLS: Epithelial cells 
% called Tumor as Epith taken

indicesTumor = intersect (indicesEpcam, indicesCD45Neg);
Tumor_XY = XY((indicesTumor),:);

% Plot scatter plot to show if the cells are present on a correct spot
scatter(-Tumor_XY(:,1), -Tumor_XY(:,2), 'filled', 'MarkerEdgeColor',[0    0.4471    0.7412],'MarkerFaceColor',[0    0.4471    0.7412] );
%scatter(-Tumor_XY(:,1), -Tumor_XY(:,2), 'filled');
hold on 
%xlim ([-3500 -1700]);
%ylim ([-3500 -1900]);
xlim ([-5000 0]);
ylim ([-4000 0]);

%% CELLS: Leukocytes
% don't match the fcs gating as no option for free shape of gate. Will be
% underepresented in Panobinostat condition mainly

indicesLeukocytes = intersect (indicesEpcamNeg, indicesCD45);
Leukocytes_XY = XY((indicesLeukocytes),:);

% Plot scatter plot to show if the cells are present on a correct spot
scatter(-Leukocytes_XY(:,1), -Leukocytes_XY(:,2), 'filled', 'MarkerEdgeColor',[0.4667    0.6745    0.1882],'MarkerFaceColor',[0.4667    0.6745    0.1882]);
hold on 
xlim ([-5000 0]);
ylim ([-4000 0]);

%% F480 Negative

F480=fcsData(:, 20);

% identify indices with values below threshold
indicesF480Neg = find(abs(F480)<fcsThreshold(18));F480(indicesF480Neg)=[];
F480Neg_XY = XY((indicesF480Neg),:); 

% Plot scatter plot to show if the cells are present on a correct spot
scatter(-F480Neg_XY(:,1), -F480Neg_XY(:,2),'filled');
hold on 
xlim ([-5000 0]);
ylim ([-5000 0]);

%% F480 positive

% delete F480 variable
F480=fcsData(:, 20);

% identify indices with values above threshold
indicesF480 = find(abs(F480)>fcsThreshold(18));F480(indicesF480)=[];
F480_XY = XY((indicesF480),:); 

% Plot scatter plot to show if the cells are present on a correct spot
scatter(-F480_XY(:,1), -F480_XY(:,2));
hold on 
xlim ([-5000 0]);
ylim ([-4000 0]);

%% CD11c negative cells

CD11c=fcsData(:, 8);

% identify indices with values below threshold
indicesCD11cNeg = find(abs(CD11c)<fcsThreshold(6));CD11c(indicesCD11cNeg)=[];
CD11cNeg_XY = XY((indicesCD11cNeg),:); 

% Plot frequency of expresion rate and scatter plot to show if the cells are present on a correct spot
scatter(-CD11cNeg_XY(:,1), -CD11cNeg_XY(:,2),'filled');
xlim ([-5000 0]);
ylim ([-4000 0]);

%% CD11c positive

CD11c=fcsData(:, 8);

% identify indices with values above threshold
indicesCD11c = find(abs(CD11c)>fcsThreshold(6));CD11c(indicesCD11c)=[];
CD11c_XY = XY((indicesCD11c),:); 

% Plot frequency of expresion rate and scatter plot to show if the cells are present on a correct spot
scatter(-CD11c_XY(:,1), -CD11c_XY(:,2),'filled');
hold on 
%xlim ([-5000 0]);
%ylim ([-4000 0]);

%% MHCII negative

MHCII=fcsData(:, 30);

% identify indices with values above threshold
indicesMHCIINeg = find(abs(MHCII)<fcsThreshold(29));MHCII(indicesMHCIINeg)=[];
MHCIINeg_XY = XY((indicesMHCIINeg),:); 

% Plot scatter plot to show if the cells are present on a correct spot
scatter(-MHCIINeg_XY(:,1), -MHCIINeg_XY(:,2));
hold on 
xlim ([-5000 0]);
ylim ([-4000 0]);
%% MHCII positive

MHCII=fcsData(:, 30);

% identify indices with values above threshold
indicesMHCII = find(abs(MHCII)>fcsThreshold(28));MHCII(indicesMHCII)=[];
MHCII_XY = XY((indicesMHCII),:); 

% Plot scatter plot to show if the cells are present on a correct spot
scatter(-MHCII_XY(:,1), -MHCII_XY(:,2));
hold on 
xlim ([-5000 0]);
ylim ([-5000 0]);

%% CSF1R positive
CSF1R=fcsData(:, 15);

% identify indices with values above threshold
indicesCSF1R = find(abs(CSF1R)>fcsThreshold(13));CSF1R(indicesCSF1R)=[];
CSF1R_XY = XY((indicesCSF1R),:); 

% Plot frequency of expresion rate and scatter plot to show if the cells are present on a correct spot
scatter(-CSF1R_XY(:,1), -CSF1R_XY(:,2),'filled');
hold on 
xlim ([-5000 0]);
ylim ([-4000 0]);

%% CELLS: Macrophages on CD45
% for display skip Epcam dependence vs FCS-based quantification

indicesMacrophages = intersect (intersect (indicesCD45, indicesF480), indicesCD11cNeg);
Macrophages_XY = XY((indicesMacrophages),:);

% Plot scatter plot to show if the cells are present on a correct spot
scatter(-Macrophages_XY(:,1), -Macrophages_XY(:,2), 'filled', 'MarkerEdgeColor',[0.0627    0.1961    0.4706],'MarkerFaceColor',[0.0627    0.1961    0.4706] );
hold on 
xlim ([-5000 0]);
ylim ([-5000 0]);

%z5=intersect(z4,intersect(intersect(z1,z2),z3))

%% CELLS: Antigen presenting Macrophages on CD45

indicesAPM = intersect (indicesMHCII, intersect (intersect (indicesCD45, indicesF480), indicesCD11cNeg));
APM_XY = XY((indicesAPM),:);

% Plot scatter plot to show if the cells are present on a correct spot
scatter(-APM_XY(:,1), -APM_XY(:,2), 'filled', 'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0]);
%scatter(-APM_XY(:,1), -APM_XY(:,2), 'filled');
hold on 
%xlim ([-3500 -1700]);
%ylim ([-3500 -1900]);
xlim ([-5000 0]);
ylim ([-4000 0]);

%% CELLS: M2 CSF1R positive Macrophages

indicesM2M= intersect (indicesMacrophages, indicesCSF1R);
M2M_XY = XY((indicesM2M),:);

% Plot scatter plot to show if the cells are present on a correct spot
scatter(-M2M_XY(:,1), -M2M_XY(:,2), 'filled', 'MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[0 0 1]);
%scatter(-M2M_XY(:,1), -M2M_XY(:,2), 'filled');
hold on 
xlim ([-5000 0]);
ylim ([-4000 0]);

%% CELLS: Dendritic cells on CD45

indicesDCs = intersect (intersect (indicesCD45, indicesCD11c), indicesF480Neg); %indicesLeukocytes instead %CD45
DCs_XY = XY((indicesDCs),:);

% Plot scatter plot to show if the cells are present on a correct spot
%scatter(-DCs_XY(:,1), -DCs_XY(:,2), 'filled', 'MarkerEdgeColor',[0.4667    0.6745    1],'MarkerFaceColor',[0.4667    0.6745    1]);
%scatter(-DCs_XY(:,1), -DCs_XY(:,2), 'filled','MarkerEdgeColor',[0.7176    0.2745    1.0000],'MarkerFaceColor',[0.7176    0.2745    1.0000] );
scatter(-DCs_XY(:,1), -DCs_XY(:,2), 'filled');
hold on 
xlim ([-5000 0]);
ylim ([-4000 0]);

%% CELLS: Antigen presenting Dendritic cells on CD45

indicesAPDCs = intersect (indicesMHCII, intersect (intersect (indicesCD45, indicesF480Neg), indicesCD11c));
APDCs_XY = XY((indicesAPDCs),:);

% Plot scatter plot to show if the cells are present on a correct spot
scatter(-APDCs_XY(:,1), -APDCs_XY(:,2), 'filled', 'MarkerEdgeColor',[0.3020    0.7451    0.9333],'MarkerFaceColor',[0.3020    0.7451    0.9333]);
hold on 
xlim ([-5000 0]);
ylim ([-5000 0]);

%% CELLS: Dendritic cells on Epithelial

indicesEpDCs = intersect (intersect (indicesTumor, indicesCD11c), indicesF480Neg); %indicesLeukocytes instead %CD45
EpDCs_XY = XY((indicesEpDCs),:);

% Plot scatter plot to show if the cells are present on a correct spot
scatter(-EpDCs_XY(:,1), -EpDCs_XY(:,2), 'filled', 'MarkerEdgeColor',[0.9294    0.6941    0.1255],'MarkerFaceColor',[0.9294    0.6941    0.1255]);
hold on 
xlim ([-5000 0]);
ylim ([-5000 0]);
%% RUN CELLS: Negative for CD11c and F480

indicesNegative = intersect (intersect (indicesCD45, indicesCD11cNeg), indicesF480Neg); 
Negative_XY = XY((indicesNegative),:);

% Plot scatter plot to show if the cells are present on a correct spot
%scatter(-Negative_XY(:,1), -Negative_XY(:,2), 'filled');
%hold on 
%xlim ([-5000 0]);
%ylim ([-4000 0]);

%% CD11b Positive
CD11b=fcsData(:, 7);

% identify indices with values above threshold
indicesCD11b = find(abs(CD11b)>fcsThreshold(5));CD11b(indicesCD11b)=[];
CD11b_XY = XY((indicesCD11b),:); 

% Plot frequency of expresion rate and scatter plot to show if the cells are present on a correct spot
scatter(-CD11b_XY(:,1), -CD11b_XY(:,2));
hold on 
xlim ([-5000 0]);
ylim ([-5000 0]);

%% Ly6G negative
Ly6G=fcsData(:, 29);

% identify indices with values above threshold
indicesLy6GNeg = find(abs(Ly6G)<fcsThreshold(27));Ly6G(indicesLy6GNeg)=[];
Ly6GNeg_XY = XY((indicesLy6GNeg),:); 

% Plot frequency of expresion rate and scatter plot to show if the cells are present on a correct spot
scatter(-Ly6GNeg_XY(:,1), -Ly6GNeg_XY(:,2));
hold on 
xlim ([-5000 0]);
ylim ([-4000 0]);

%% Ly6G positive
Ly6G=fcsData(:, 29);

% identify indices with values above threshold
indicesLy6G = find(abs(Ly6G)>fcsThreshold(27));Ly6G(indicesLy6G)=[];
Ly6G_XY = XY((indicesLy6G),:); 

% Plot frequency of expresion rate and scatter plot to show if the cells are present on a correct spot
scatter(-Ly6G_XY(:,1), -Ly6G_XY(:,2));
hold on 
xlim ([-5000 0]);
ylim ([-5000 0]);

%% CELLS: Neutrophils

indicesNeutrophils = intersect (intersect (indicesNegative, indicesCD11b), indicesLy6G); 
Neutrophils_XY = XY((indicesNeutrophils),:);

% Plot scatter plot to show if the cells are present on a correct spot
scatter(-Neutrophils_XY(:,1), -Neutrophils_XY(:,2), 'filled', 'MarkerEdgeColor',[0.9294    0.6941    0.1255],'MarkerFaceColor',[0.9294    0.6941    0.1255]);
%scatter(-Neutrophils_XY(:,1), -Neutrophils_XY(:,2), 'filled');
hold on 
%xlim ([-3500 -1700]);
%ylim ([-3500 -1900]);
xlim ([-5000 0]);
ylim ([-4000 0]);

%% CELLS: MHCII positive Neutrophils

indicesAPN= intersect (indicesNeutrophils, indicesMHCII);
APN_XY = XY((indicesAPN),:);

% Plot scatter plot to show if the cells are present on a correct spot
scatter(-APN_XY(:,1), -APN_XY(:,2), 'filled', 'MarkerEdgeColor',[0.0824    0.4000    0.1059],'MarkerFaceColor',[0.0824    0.4000    0.1059]);
%scatter(-APN_XY(:,1), -APN_XY(:,2), 'filled');
hold on 
%xlim ([-3500 -1700]);
%ylim ([-3500 -1900]);
xlim ([-5000 0]);
ylim ([-4000 0]);

%% CELLS: Immature myeloid cells

indicesImmature = intersect (intersect (indicesNegative, indicesCD11b), indicesLy6GNeg); 
Immature_XY = XY((indicesImmature),:);

% Plot scatter plot to show if the cells are present on a correct spot
scatter(-Immature_XY(:,1), -Immature_XY(:,2), 'filled', 'MarkerEdgeColor',[0.5020    0.5020    0.5020],'MarkerFaceColor',[0.5020    0.5020    0.5020]);
scatter(-Immature_XY(:,1), -Immature_XY(:,2), 'filled');
hold on 
%xlim ([-3500 -1700]);
%ylim ([-3500 -1900]);
xlim ([-5000 0]);
ylim ([-4000 0]);

%% CD31 Negative
CD31=fcsData(:, 9);

% identify indices with values above threshold
indicesCD31Neg = find(abs(CD31)<fcsThreshold(7));CD31(indicesCD31Neg)=[];
CD31Neg_XY = XY((indicesCD31Neg),:); 

% Plot frequency of expresion rate and scatter plot to show if the cells are present on a correct spot
scatter(-CD31Neg_XY(:,1), -CD31Neg_XY(:,2));
hold on 
xlim ([-5000 0]);
ylim ([-4000 0]);

%% CD31 Positive
CD31=fcsData(:, 9);

% identify indices with values above threshold
indicesCD31 = find(abs(CD31)>fcsThreshold(7));CD31(indicesCD31)=[];
CD31_XY = XY((indicesCD31),:); 

% Plot frequency of expresion rate and scatter plot to show if the cells are present on a correct spot
scatter(-CD31_XY(:,1), -CD31_XY(:,2));
hold on 
xlim ([-5000 0]);
ylim ([-4000 0]);

%% aSMA Negative
aSMA=fcsData(:, 4);

% identify indices with values above threshold
indicesaSMANeg = find(abs(aSMA)<fcsThreshold(2));aSMA(indicesaSMANeg)=[];
aSMANeg_XY = XY((indicesaSMANeg),:); 

% Plot frequency of expresion rate and scatter plot to show if the cells are present on a correct spot
scatter(-aSMANeg_XY(:,1), -aSMANeg_XY(:,2), 'filled');
hold on 
xlim ([-5000 0]);
ylim ([-4000 0]);

%% aSMA Positive
aSMA=fcsData(:, 4);

% identify indices with values above threshold
indicesaSMA = find(abs(aSMA)>fcsThreshold(2));aSMA(indicesaSMA)=[];
aSMA_XY = XY((indicesaSMA),:); 

% Plot frequency of expresion rate and scatter plot to show if the cells are present on a correct spot
scatter(-aSMA_XY(:,1), -aSMA_XY(:,2), 'filled');
hold on 
xlim ([-5000 0]);
ylim ([-5000 0]);

%% CELLS: Endothelial Cells

indicesEndothelial = intersect (indicesCD31, intersect (intersect (indicesCD45Neg, indicesEpcamNeg), indicesaSMANeg));
Endothelial_XY = XY((indicesEndothelial),:);

% Plot scatter plot to show if the cells are present on a correct spot
%scatter(-Endothelial_XY(:,1), -Endothelial_XY(:,2), 'filled', 'MarkerEdgeColor',[0.6353    0.0784    0.1843],'MarkerFaceColor',[0.6353    0.0784    0.1843]);
scatter(-Endothelial_XY(:,1), -Endothelial_XY(:,2), 'filled');
hold on 
xlim ([-5000 0]);
ylim ([-4000 0]);

%% CELLS: Mesenchymal Cells

indicesMesenchymal = intersect (indicesaSMA, intersect (intersect (indicesCD45Neg, indicesEpcamNeg), indicesCD31Neg));
Mesenchymal_XY = XY((indicesMesenchymal),:);

% Plot scatter plot to show if the cells are present on a correct spot
scatter(-Mesenchymal_XY(:,1), -Mesenchymal_XY(:,2), 'filled', 'MarkerEdgeColor',[0.4941    0.1843    0.5569],'MarkerFaceColor',[0.4941    0.1843    0.5569]);
%scatter(-Mesenchymal_XY(:,1), -Mesenchymal_XY(:,2), 'filled');
hold on 
%xlim ([-3500 -1700]);
%ylim ([-3500 -1900]);
xlim ([-5000 0]);
ylim ([-4000 0]);

%% CELLS: Pericytes

indicesPericytes = intersect (indicesaSMA, intersect (intersect (indicesCD45Neg, indicesEpcamNeg), indicesCD31));
pericytes_XY = XY((indicesPericytes),:);

% Plot scatter plot to show if the cells are present on a correct spot
scatter(-pericytes_XY(:,1), -pericytes_XY(:,2), 'filled', 'MarkerEdgeColor',[0.8510    0.3255    0.0980],'MarkerFaceColor',[0.8510    0.3255    0.0980]);
%scatter(-CSC_XY(:,1), -CSC_XY(:,2), 'filled');
hold on 
xlim ([-5000 0]);
ylim ([-4000 0]);


%% CELLS: non-immune stroma

indicesNonIm = union (indicesPericytes, (union (indicesMesenchymal, indicesEndothelial)));
NonIm_XY = XY((indicesNonIm),:);

% Plot scatter plot to show if the cells are present on a correct spot
scatter(-NonIm_XY(:,1), -NonIm_XY(:,2), 'filled', 'MarkerEdgeColor',[0.8510    0.3255    0.0980],'MarkerFaceColor',[0.8510    0.3255    0.0980]);
%scatter(-CSC_XY(:,1), -CSC_XY(:,2), 'filled');
hold on 
xlim ([-5000 0]);
ylim ([-4000 0]);

%% CELLS: EndothPericytes

indicesEndPer = union (indicesPericytes, indicesEndothelial);
EndPer_XY = XY((indicesEndPer),:);

% Plot scatter plot to show if the cells are present on a correct spot
%scatter(-EndPer_XY(:,1), -EndPer_XY(:,2), 'filled', 'MarkerEdgeColor',[0.8510    0.3255    0.0980],'MarkerFaceColor',[0.8510    0.3255    0.0980]);
scatter(-EndPer_XY(:,1), -EndPer_XY(:,2), 'filled');
hold on 
xlim ([-5000 0]);
ylim ([-4000 0]);
%% Sox9
Sox9=fcsData(:, 37);

% identify indices with values above threshold
indicesSox9 = find(abs(Sox9)>fcsThreshold(35));Sox9(indicesSox9)=[];
Sox9_XY = XY((indicesSox9),:); 

% Plot frequency of expresion rate and scatter plot to show if the cells are present on a correct spot
scatter(-Sox9_XY(:,1), -Sox9_XY(:,2));
hold on 
xlim ([-5000 0]);
ylim ([-5000 0]);

%% Ki67 Negative
Ki67=fcsData(:, 27);

% identify indices with values above threshold
indicesKi67Neg = find(abs(Ki67)<fcsThreshold(25));Ki67(indicesKi67Neg)=[];
Ki67Neg_XY = XY((indicesKi67Neg),:); 

% Plot frequency of expresion rate and scatter plot to show if the cells are present on a correct spot
scatter(-Ki67Neg_XY(:,1), -Ki67Neg_XY(:,2));
hold on 
xlim ([-5000 0]);
ylim ([-4000 0]);

%% Ki67 positive
Ki67=fcsData(:, 27);

% identify indices with values above threshold
indicesKi67 = find(abs(Ki67)>fcsThreshold(25));Ki67(indicesKi67)=[];
Ki67_XY = XY((indicesKi67),:); 

% Plot frequency of expresion rate and scatter plot to show if the cells are present on a correct spot
scatter(-Ki67_XY(:,1), -Ki67_XY(:,2));
hold on 
xlim ([-5000 0]);
ylim ([-4000 0]);

%% PyMT
PyMT=fcsData(:, 19);
% identifies tumor (vs all epithelial) cells in the MMTV-PyMT model

% identify indices with values above threshold
indicesPyMT = find(abs(PyMT)>fcsThreshold(17));PyMT(indicesPyMT)=[];
PyMT_XY = XY((indicesPyMT),:); 

% Plot frequency of expresion rate and scatter plot to show if the cells are present on a correct spot
scatter(-PyMT_XY(:,1), -PyMT_XY(:,2));
hold on 
xlim ([-5000 0]);
ylim ([-4000 0]);

%% CELLS: Cancer Stem cells

indicesCSC = intersect (indicesSox9, intersect (intersect (indicesTumor, indicesKi67Neg), indicesPyMT));
CSC_XY = XY((indicesCSC),:);

% Plot scatter plot to show if the cells are present on a correct spot
%scatter(-CSC_XY(:,1), -CSC_XY(:,2), 'filled', 'MarkerEdgeColor',[0.4667    0.6745    0.1882],'MarkerFaceColor',[0.4667    0.6745    0.1882]);
scatter(-CSC_XY(:,1), -CSC_XY(:,2), 'filled');
hold on 
%xlim ([-3500 -1700]);
%ylim ([-3500 -1900]);
xlim ([-5000 0]);
ylim ([-4000 0]);

%% CD3
CD3=fcsData(:, 10);

% identify indices with values above threshold
indicesCD3 = find(abs(CD3)>fcsThreshold(8));CD3(indicesCD3)=[];
CD3_XY = XY((indicesCD3),:); 

% Plot frequency of expresion rate and scatter plot to show if the cells are present on a correct spot

scatter(-CD3_XY(:,1), -CD3_XY(:,2));
hold on 
title ("CD3 position in ROI", 'Fontsize', 14);
hold on 
xlim ([-5000 0]);
ylim ([-4000 0]);
%% CELLS:  T cells

indicesCD3T = intersect (intersect (indicesCD45, indicesCD3), indicesEpcamNeg);
CD3T_XY = XY((indicesCD3T),:);

% Plot scatter plot to show if the cells are present on a correct spot
scatter(-CD3T_XY(:,1), -CD3T_XY(:,2), 'filled', 'MarkerEdgeColor',[1.0000    0.0745    0.6510],'MarkerFaceColor',[1.0000    0.0745    0.6510] );
hold on 
xlim ([-5000 0]);
ylim ([-4000 0]);

%% CD4
CD4=fcsData(:, 13);

% identify indices with values above threshold
indicesCD4 = find(abs(CD4)>fcsThreshold(11));CD4(indicesCD4)=[];
CD4_XY = XY((indicesCD4),:); 

% Plot frequency of expresion rate and scatter plot to show if the cells are present on a correct spot
scatter(-CD4_XY(:,1), -CD4_XY(:,2), 'filled');
hold on 
title ("CD4 position in ROI", 'Fontsize', 14);
hold on 
xlim ([-5000 0]);
ylim ([-4000 0]);

%% CELLS: CD4 T cells

% not many cells intratumorally in the MIMA study. Even if Epcam, CD8 independent 

indicesCD4T = intersect (intersect (indicesCD45, indicesCD3T), indicesCD4);
CD4T_XY = XY((indicesCD4T),:);

% Plot scatter plot to show if the cells are present on a correct spot
scatter(-CD4T_XY(:,1), -CD4T_XY(:,2), 'filled', 'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0] );
hold on 
xlim ([-5000 0]);
ylim ([-4000 0]);

%% CD8
CD8=fcsData(:, 14);

% identify indices with values above threshold
indicesCD8 = find(abs(CD8)>fcsThreshold(12));CD8(indices)=[];
CD8_XY = XY((indicesCD8),:); 
% Plot frequency of expresion rate and scatter plot to show if the cells are present on a correct spot

scatter(-CD8_XY(:,1), -CD8_XY(:,2), 'filled');
hold on 
title ("CD8 position in ROI", 'Fontsize', 14);
hold on 
xlim ([-4000 0]);
ylim ([-4000 0]);

%% CELLS: CD8 T cells

indicesCD8T = intersect (intersect (indicesCD45, indicesCD3T), indicesCD8);
CD8T_XY = XY((indicesCD8T),:);

% Plot scatter plot to show if the cells are present on a correct spot
scatter(-CD8T_XY(:,1), -CD8T_XY(:,2), 'filled', 'MarkerEdgeColor',[0.0588    1.0000    1.0000],'MarkerFaceColor',[0.0588    1.0000    1.0000] );
hold on 
xlim ([-5000 0]);
ylim ([-4000 0]);

%% CELLS: Negative for all

indicesAllNeg = intersect (indicesEpcamNeg, intersect (intersect (indicesCD45Neg, indicesCD31Neg), indicesaSMANeg));
AllNeg_XY = XY((indicesAllNeg),:);

% Plot scatter plot to show if the cells are present on a correct spot
scatter(-AllNeg_XY(:,1), -AllNeg_XY(:,2), 'filled', 'MarkerEdgeColor',[0.4667    0.6745    0.1882],'MarkerFaceColor',[0.4667    0.6745    0.1882]);
hold on 
%xlim ([-3500 -1700]);
%ylim ([-3500 -1900]);
xlim ([-5000 0]);
ylim ([-4000 0]);

%% CELLS: Proliferating Tumor

indicesTKi67 = intersect (intersect (indicesTumor, indicesKi67), indicesPyMT);
TKi67_XY = XY((indicesTKi67),:);

% Plot scatter plot to show if the cells are present on a correct spot
scatter(-TKi67_XY(:,1), -TKi67_XY(:,2), 'filled', 'MarkerEdgeColor',[0.3020    0.7451    0.9333],'MarkerFaceColor',[0.3020    0.7451    0.9333]);
%scatter(-TKi67_XY(:,1), -TKi67_XY(:,2), 'filled');
hold on 
%xlim ([-3500 -1700]);
%ylim ([-3500 -1900]);
xlim ([-5000 0]);
ylim ([-4000 0]);

%% find any cell type of interest for display and quatification

%% CALR
CALR=fcsData(:, 5);

% identify indices with values above threshold
indicesCALR = find(abs(CALR)>fcsThreshold(3));CALR(indicesCALR)=[];
CALR_XY = XY((indicesCALR),:); 
% Plot frequency of expresion rate and scatter plot to show if the cells are present on a correct spot

scatter(-CALR_XY(:,1), -CALR_XY(:,2), 'filled');
hold on 
title ("CALR Sox9 position in ROI", 'Fontsize', 14);
hold on 
xlim ([-5000 0]);
ylim ([-4000 0]);

%% Gal3
Gal3=fcsData(:, 22);

% identify indices with values above threshold
indices = find(abs(Gal3)>fcsThreshold(20));Gal3(indices)=[];
Gal3_XY = XY((indices),:); 
% Plot frequency of expresion rate and scatter plot to show if the cells are present on a correct spot


scatter(-Gal3_XY(:,1), -Gal3_XY(:,2), 'filled');
hold on 
title ("Galectin-3 position in ROI", 'Fontsize', 14);

%% ICAM1
ICAM1=fcsData(:, 24);

% identify indices with values above threshold
indicesICAM1 = find(abs(ICAM1)>fcsThreshold(22));ICAM1(indicesICAM1)=[];
ICAM1_XY = XY((indicesICAM1),:); 
% Plot frequency of expresion rate and scatter plot to show if the cells are present on a correct spot

scatter(-ICAM1_XY(:,1), -ICAM1_XY(:,2),'filled');
hold on 
title ("ICAM-1 position in ROI", 'Fontsize', 14);

hold on 
xlim ([-5000 0]);
ylim ([-4000 0]);
%% CELLS: ICAM1 on tumor cells

indicesICAM1_T = intersect (intersect (indicesTumor, indicesPyMT), indicesICAM1); %indicesLeukocytes instead %CD45
ICAM1_T_XY = XY((indicesICAM1_T),:);

% Plot scatter plot to show if the cells are present on a correct spot
%scatter(-DCs_XY(:,1), -DCs_XY(:,2), 'filled', 'MarkerEdgeColor',[0.4667    0.6745    1],'MarkerFaceColor',[0.4667    0.6745    1]);
scatter(-ICAM1_T_XY(:,1), -ICAM1_T_XY(:,2), 'filled');
hold on 
xlim ([-5000 0]);
ylim ([-4000 0]);

%% MHCI
MHCI=fcsData(:, 31);

% identify indices with values above threshold
indicesMHCI = find(abs(MHCI)>fcsThreshold(29));MHCI(indicesMHCI)=[];
MHCI_XY = XY((indicesMHCI),:); 
% Plot frequency of expresion rate and scatter plot to show if the cells are present on a correct spot

%scatter(-MHCI_XY(:,1), -MHCI_XY(:,2),'filled');
%hold on 
%title ("MHC-I high position in ROI", 'Fontsize', 14);

%% CELLS: MHC-I on tumor cells

indicesMHCI_T = intersect (intersect (indicesTumor, indicesPyMT), indicesMHCI); %indicesLeukocytes instead %CD45
MHCI_T_XY = XY((indicesMHCI_T),:);

% Plot scatter plot to show if the cells are present on a correct spot
%scatter(-DCs_XY(:,1), -DCs_XY(:,2), 'filled', 'MarkerEdgeColor',[0.4667    0.6745    1],'MarkerFaceColor',[0.4667    0.6745    1]);
scatter(-MHCI_T_XY(:,1), -MHCI_T_XY(:,2), 'filled');
hold on 
xlim ([-5000 0]);
ylim ([-4000 0]);


%% CELLS: Calr on stem cells

% Q: are CSC in CSC/DC niche high in calreticulin?

indicesCalrCSC = intersect (indicesCSC, indicesCALR); %indicesLeukocytes instead %CD45
CalrCSC_XY = XY((indicesCalrCSC),:);

% Plot scatter plot to show if the cells are present on a correct spot
%scatter(-DCs_XY(:,1), -DCs_XY(:,2), 'filled', 'MarkerEdgeColor',[0.4667    0.6745    1],'MarkerFaceColor',[0.4667    0.6745    1]);
scatter(-CalrCSC_XY(:,1), -CalrCSC_XY(:,2), 'filled');
hold on 
xlim ([-5000 0]);
ylim ([-4000 0]);


%% CELLS: N1 neutrophils

indicesN1= intersect (intersect (indicesNeutrophils, indicesArg1Neg), indicesICAM1);
N1_XY = XY((indicesN1),:);

% Plot scatter plot to show if the cells are present on a correct spot
%scatter(-APN_XY(:,1), -APN_XY(:,2), 'filled', 'MarkerEdgeColor',[0.3020    0.7451    0.9333],'MarkerFaceColor',[0.3020    0.7451    0.9333]);
scatter(-N1_XY(:,1), -N1_XY(:,2), 'filled');
hold on 
%xlim ([-3500 -1700]);
%ylim ([-3500 -1900]);
xlim ([-5000 0]);
ylim ([-4000 0]);