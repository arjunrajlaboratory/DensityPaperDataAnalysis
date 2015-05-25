%Summed dapi intensity

contents = dir('data*');
m = 1;
color = 'tmr';
clear intensityData;

for i = 1:numel(contents)
    load(contents(i).name);
    for j = 1:numel(objects)
        
        obj = objects(j);
        
        dapiMask = obj.channels.dapi.processor.mask;
        
        tmrStk = obj.channelStk('tmr');
        tmrStk = max(tmrStk,[],3);
        tmrMask = maskWithDapi(tmrStk);
        %tmrMask = bwareaopen(tmrMask1,2000); 
        
        mask = obj.object_mask.mask;

        alexaStk = obj.channelStk(color);
        
        maxAlexa = max(alexaStk,[],3);
        maxAlexa = double(maxAlexa);
        
        area = numel(find(tmrMask==1));
        
        int = maxAlexa(tmrMask == 1);
        bkg = maxAlexa(tmrMask == 0 & dapiMask == 1);
        
        nucIntens = sum(int(:));
        bkgrd = area*median(bkg(:));
        
        avgNucIntens = mean(int(:));
        avgBkgrd = median(bkg(:));
        avgNucIntens = avgNucIntens - avgBkgrd;
        
        nucIntens = nucIntens - bkgrd;
        
        intensityData(m,:) = [i j obj.channels.nir.numSpots avgNucIntens];
        m = m+1;
    end
end

dlmwrite('rRNAprecursorIntensity.txt',intensityData,'\t');