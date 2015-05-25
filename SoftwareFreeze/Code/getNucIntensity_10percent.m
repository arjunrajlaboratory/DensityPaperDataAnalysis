function [sumIntensity avgIntensity] = getNucIntensity(obj,nucColor,gapdhColor)

% Usage: [sumIntensity avgIntensity] = getNucIntensity(obj,nucColor,gapdhColor)
% Returns the mean and average intensity of the nucleus (color of your
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
crop=crop.*uint16(obj.object_mask.mask);

bkgrd = min(crop(find(crop(:))));

crop = crop-bkgrd;

nucCrop = crop.*uint16(obj.channels.dapi.processor.mask);

nonzeroCrop = nucCrop(nucCrop>0);
cutoff = quantile(nonzeroCrop,.9);
finalCrop = nonzeroCrop(nonzeroCrop>cutoff);

sumIntensity = sum(finalCrop(:));
avgIntensity = sumIntensity/length(finalCrop);