function [sumIntensity avgIntensity] = getCellIntensity(obj,nucColor,gapdhColor)

% Usage: [sumIntensity avgIntensity] = getCellIntensity(obj,nucColor,gapdhColor)
% Returns the mean and average intensity of the cell, nucleus excluded (color of your
% choice).  

minPlane = min(obj.channels.(gapdhColor).spotCoordinates(:,3));
maxPlane = max(obj.channels.(gapdhColor).spotCoordinates(:,3));
cropPlane = round((maxPlane-minPlane)/2);

if cropPlane == 0
    sumIntensity = 0;
    avgIntensity = 0;
    return;
end

img=obj.channelStk(nucColor);
crop=img(:,:,cropPlane);
bkgrd = min(crop(find(crop(:))));

crop=crop.*uint16(obj.object_mask.mask);
crop = crop-bkgrd;
crop = crop.*uint16(~obj.channels.dapi.processor.mask);

sumIntensity = sum(crop(:));
avgIntensity = sumIntensity/(numel(find(obj.object_mask.mask))-numel(find(obj.channels.dapi.processor.mask)));