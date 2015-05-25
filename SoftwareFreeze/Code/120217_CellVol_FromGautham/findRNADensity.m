function [cellVol numCountRNA numDenseRNA density celltop cellbottom] = findRNADensity( obj, denseChannel, countChannel );
%Finds density of spots in a given cell.

%Find DAPI mask to exclude volume/spots in the nucleus
dapiStk = obj.channelStk('dapi');
dapiStk = max(dapiStk,[],3);
dapiMask1 = maskWithDapi(dapiStk);
dapiMask = bwareaopen(dapiMask1,2000);

[cellVol cellHeight celltop cellbottom] = objectvolume_fromRNAspots_Olivia(obj,denseChannel,dapiMask);

numCountRNA = numel(obj.channels.(countChannel).fitdataRNAonly.xp_fit);
numDenseRNA = numel(obj.channels.(denseChannel).fitdataRNAonly.xp_fit);

x = round(obj.channels.(countChannel).fitdataRNAonly.xp_fit);
y = round(obj.channels.(countChannel).fitdataRNAonly.yp_fit);

xD = round(obj.channels.(denseChannel).fitdataRNAonly.xp_fit);
yD = round(obj.channels.(denseChannel).fitdataRNAonly.yp_fit);

ind = sub2ind(size(dapiMask),x,y);
numExclude = numel(find(dapiMask(ind)==1));

numCountRNA = numCountRNA - numExclude;

ind = sub2ind(size(dapiMask),xD,yD);
numExclude = numel(find(dapiMask(ind)==1));

numDenseRNA = numDenseRNA - numExclude;


density = numCountRNA/cellVol;

%figure; plotFISHSpots(obj, 30, 'DAPIFlag',true);

%Plot mask with dapi overlay to check
mask = obj.object_mask.mask;
mask(find(dapiMask==1))=0;
%figure; imagesc(mask);

end