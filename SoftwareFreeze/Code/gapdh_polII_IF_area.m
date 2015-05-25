% GAPDH outside of nucleus (nir)
% total cell area
% Pol II (tmr)
% Pol II IF intensity (cy)

clear data;
contents = dir('data*');
fillChannel = 'nir';
ifChannel = 'cy';
countChannel = 'tmr';

m=1;

for i = 1:numel(contents)
    load(contents(i).name);
    for j = 1:numel(objects)
        
        obj = objects(j);
        
        mask = obj.object_mask.mask;
        dapiMask = obj.channels.dapi.processor.mask;
        
        %number GAPDH outside of nucleus
        numFillRNA = size(obj.channels.(fillChannel).spotCoordinates,1);
        inNuc = numel(find(obj.channels.(fillChannel).inNucleus(dapiMask)==1));
        numFillRNA = numFillRNA - inNuc;
        
        %total cell area
        area = numel(find(mask));
        
        %Pol II mRNA
        numCountRNA = size(obj.channels.(countChannel).spotCoordinates,1);
        
        %Pol II IF
        img=obj.channelStk(ifChannel);
        cropPlane = 15;
        crop=img(:,:,cropPlane);
        crop=crop.*uint16(mask);
        
        bkgrd = min(crop(find(crop(:))));
        
        crop = crop-bkgrd;
        
        nucCrop = crop.*uint16(obj.channels.dapi.processor.mask);
        
        dapiSum = sum(nucCrop(:));
        dapiAvg = dapiSum/numel(find(obj.channels.dapi.processor.mask));
        
        data(m,:) = [i j area numFillRNA numCountRNA dapiSum dapiAvg];
        m = m+1;
    end
end

dlmwrite('gapdh_polII_IF_area.txt',data,'\t');