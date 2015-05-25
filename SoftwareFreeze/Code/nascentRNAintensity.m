contents = dir('data*');
m = 1;
color = 'alexa';
clear intensityData;

for i = 1:numel(contents)
    load(contents(i).name);
    for j = 1:numel(objects)
        
        %if i==4 & j==1
        %    continue;
        %end
        
        obj = objects(j);
        
        dapiMask = obj.channels.dapi.processor.mask;
        mask = obj.object_mask.mask;

        alexaStk = obj.channelStk(color);
        
        maxAlexa = max(alexaStk,[],3);
        maxAlexa = double(maxAlexa);
        
        area = numel(find(dapiMask==1));
        
        int = maxAlexa(dapiMask == 1);
        bkg = maxAlexa(dapiMask == 0 & mask == 1);
        
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

dlmwrite('alexaIntensity.txt',intensityData,'\t');