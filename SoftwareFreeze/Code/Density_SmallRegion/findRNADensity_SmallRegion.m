function [cellvolume numCountRNA numDenseRNA density celltop cellbottom] = findRNADensity_SmallRegion( obj, denseChannel, countChannel );
%Finds density of spots in a given cell.

%Find DAPI mask to exclude volume/spots in the nucleus
dapiStk = obj.channelStk('dapi');
dapiStk = max(dapiStk,[],3);
dapiMask1 = maskWithDapi(dapiStk);
dapiMask = bwareaopen(dapiMask1,2000);

%Use this to find new convex mask inside of cell
mask = obj.object_mask.mask;
mask(find(dapiMask==1))=0;
figure; imagesc(mask);

newROI = imfreehand;

newMask = newROI.createMask;

% we always have to do this x<->y shift for things to work correctly. 
xp=obj.channels.(denseChannel).fitdataRNAonly.yp_fit;
yp=obj.channels.(denseChannel).fitdataRNAonly.xp_fit;
zp=obj.channels.(denseChannel).fitdataRNAonly.rawzp;

rnacoords=[xp', yp', zp'];

%tentconstruct_fromRNAspots_forceOutline_Olivia forces "spots" on the mask
%boundary to be included in the Delaunay triangulation to get a better
%volume estimate.

[cellvolume cellheight celltop cellbottom]= tentconstruct_fromRNAspots_Olivia( dapiMask, newMask, rnacoords, obj.channels.(denseChannel).maxlaplaceimage );
%[cellvolume cellheight celltop cellbottom]= tentconstruct_fromRNAspots_forceOutline_Olivia( dapiMask, newMask , rnacoords, obj.channels.(denseChannel).maxlaplaceimage );

numCountRNA = numel(obj.channels.(countChannel).fitdataRNAonly.xp_fit);
numDenseRNA = numel(obj.channels.(denseChannel).fitdataRNAonly.xp_fit);

x = round(obj.channels.(countChannel).fitdataRNAonly.xp_fit);
y = round(obj.channels.(countChannel).fitdataRNAonly.yp_fit);

xD = round(obj.channels.(denseChannel).fitdataRNAonly.xp_fit);
yD = round(obj.channels.(denseChannel).fitdataRNAonly.yp_fit);

ind = sub2ind(size(newMask),x,y);
numExclude = numel(find(newMask(ind)==0));

numCountRNA = numCountRNA - numExclude;

ind = sub2ind(size(newMask),xD,yD);
numExclude = numel(find(newMask(ind)==0));

numDenseRNA = numDenseRNA - numExclude;


density = numCountRNA/cellvolume;

end