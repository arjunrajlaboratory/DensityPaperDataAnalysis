%Summed dapi intensity

contents = dir('data*');
m = 1;
clear intensityData;

for i = 1:numel(contents)
    load(contents(i).name);
    for j = 1:numel(objects)
        
        obj = objects(j);
        
        dapiMask = obj.channels.dapi.processor.mask;

        dapiStk = obj.channelStk('dapi');
        
        maxDapi = max(dapiStk,[],3);
        maxDapi = double(maxDapi);
        
        area = numel(find(dapiMask==1));
        
        int = maxDapi(dapiMask == 1);
        bkg = maxDapi(dapiMask == 0);
        
        nucIntens = sum(int(:));
        bkgrd = area*median(bkg(:));
        
        avgNucIntens = mean(int(:));
        avgBkgrd = median(bkg(:));
        avgNucIntens = avgNucIntens - avgBkgrd;
        
        nucIntens = nucIntens - bkgrd;
        
        intensityData(m,:) = [i j nucIntens avgNucIntens];
        m = m+1;
    end
end

dlmwrite('nucIntensity.txt',intensityData,'\t');