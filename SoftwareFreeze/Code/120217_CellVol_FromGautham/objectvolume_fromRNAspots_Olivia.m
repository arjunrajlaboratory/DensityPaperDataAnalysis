function [ cellvolume cellheight celltop cellbottom ] = objectvolume_fromRNAspots_Olivia( obj, rna_color, dapiMask )
% obj is the image object representing for example a cell. 
% rna_color is a string indicating which color the RNA spots are in.
% the image object must have had its RNA counted by the gaussian method.
% The spots in obj.channels.(rna_color).fitdataRNAonly will be used to
% generate the volume

% we always have to do this x<->y shift for things to work correctly. 
xp=obj.channels.(rna_color).fitdataRNAonly.yp_fit;
yp=obj.channels.(rna_color).fitdataRNAonly.xp_fit;
zp=obj.channels.(rna_color).fitdataRNAonly.rawzp;

rnacoords=[xp', yp', zp'];

%tentconstruct_fromRNAspots_forceOutline_Olivia forces "spots" on the mask
%boundary to be included in the Delaunay triangulation to get a better
%volume estimate.

[cellvolume cellheight celltop cellbottom]= tentconstruct_fromRNAspots_Olivia( dapiMask, obj.object_mask.mask , rnacoords, obj.channels.(rna_color).maxlaplaceimage );
%[cellvolume cellheight celltop cellbottom]= tentconstruct_fromRNAspots_forceOutline_Olivia( dapiMask, obj.object_mask.mask , rnacoords, obj.channels.(rna_color).maxlaplaceimage );

%To compare cell height maps from both methods:
%[cellvolume cellheight]= tentconstruct_fromRNAspots_Olivia( dapiMask, obj.object_mask.mask , rnacoords, obj.channels.(rna_color).maxlaplaceimage );
%[cellvolume_force cellheight_force]= tentconstruct_fromRNAspots_forceOutline_Olivia( dapiMask, obj.object_mask.mask , rnacoords, obj.channels.(rna_color).maxlaplaceimage );

%figure; imagesc(cellheight); colormap jet;
%figure; imagesc(cellheight_force); colormap jet;

end

