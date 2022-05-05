% standard cell type presentation (rate and significance) in form of a
% heatmap

colormap('pink');
imagesc (sss_msa);
%caxis ([0 0.61])
colorbar;

set(gca,'Xtick',1:16,'XTickLabel',{'leukocytes','T cells', 'CD4 T cells', 'CD8 T cells', 'Macrophages (M)', 'Protumorigenic M', 'Antigen presenting M', 'Dendtritic cells (DC)', 'Antigen presenting DC', 'Neutrophils (N)', 'Antigen presenting N', 'Immature myeloid cells', 'Endothelial cells', 'Mesenchymal cells', 'Pericytes', 'Epcam leukocytes'})
%xticklabels ({'Immediate','Proximal', 'Proximal 2', 'Proximal 3', 'Proximal 4', 'Proximal 5', 'Border', 'Border 2', 'Border 3', 'Border 4', 'Border 5', 'Border 6', 'Border 7', 'Border 8', 'Distal', 'Distal 2', 'Distal 3', 'Distal 4', 'Distal 5', 'Distal 6', 'Remote', 'Remote 2', 'Remote 3', 'Remote 4','Remote 5', 'Remote 6'});
xticklabels ({'proximal','border', 'distal', 'random control', 'Macrophages (M)', 'Protumorigenic M', 'Antigen presenting M', 'Dendtritic cells (DC)', 'Antigen presenting DC', 'Neutrophils (N)', 'Antigen presenting N', 'Immature myeloid cells', 'Endothelial cells', 'Mesenchymal cells', 'Pericytes', 'Epcam leukocytes'});
set (gca,'XTicklabelRotation', 35, 'Fontsize', 16)
%yticklabels ({'Olaparib','Lenvatinib', 'Palbociclib', 'Venetoclax', 'Panobinostat', 'Paclitaxel', 'Doxorubicin','Pano/Veneto (d8)', 'PEG control'});
set (gca, 'Fontsize', 16) %'Pano/Veneto', 
set (gca, 'Fontsize', 16)

%hold on
title ('Rate of significance', 'Fontsize', 18)



