clear;
contents = dir('data*');
fillChannel = 'nir';
countChannel = 'alexa';
%countChannel = 'tmr';
cyclinChannel = 'cy';
intChannel = 'tmr';

m=1;

for i = 1:numel(contents)
    load(contents(i).name);
    for j = 1:numel(objects)
        
        %if i==10 & j==1
        %    continue;
        %end
        
        obj = objects(j);
        mask = obj.object_mask.mask;
        
        %Find dapi mask
        %dapiStk = obj.channelStk('dapi');
        %dapiStk = max(dapiStk,[],3);
        %dapiMask1 = maskWithDapi(dapiStk);
        %dapiMask = bwareaopen(dapiMask1,2000);
        dapiMask = findDapiMask(obj);
        
        %area = numel(find((mask & ~dapiMask)==1));
        
        %areaData(m,:) = [i j area];
        
        %numFillRNA = numel(obj.metadata.(fillChannel).normal.x);
        %numCountRNA = numel(obj.metadata.(countChannel).normal.x);
        numFillRNA = numel(obj.metadata.(fillChannel).y);
        numCountRNA = numel(obj.channels.(countChannel).fitdataRNAonly.xp_fit);
        numCyclinRNA = numel(obj.channels.(cyclinChannel).fitdataRNAonly.xp_fit);
        if obj.channels.(cyclinChannel).has_clearthreshold==0
            numCyclinRNA = 0;
        end
        numIntRNA = numel(obj.channels.(intChannel).fitdataRNAonly.xp_fit);
        if obj.channels.(intChannel).has_clearthreshold==0
            numIntRNA = 0;
        end
        
        x = round(obj.channels.(countChannel).fitdataRNAonly.yp_fit);
        y = round(obj.channels.(countChannel).fitdataRNAonly.xp_fit);
        %x(find(x==0))=1;
        %y(find(y==0))=1;
        %x(find(x>size(dapiMask,2)))=size(dapiMask,2);
        %y(find(y>size(dapiMask,1)))=size(dapiMask,1);
        
        xD = round(obj.metadata.(fillChannel).y);
        yD = round(obj.metadata.(fillChannel).x);
        
        %ind = sub2ind(size(dapiMask),y,x);
        %numExclude = numel(find(dapiMask(ind)==1));
        numExclude = findSpotsInNuc(x,y,dapiMask);
        
        numCountRNA = numCountRNA - numExclude;
        
        %ind = sub2ind(size(dapiMask),xD,yD);
        %numExclude = numel(find(dapiMask(ind)==1));
        numExclude = findSpotsInNuc(yD,xD,dapiMask);
        
        numFillRNA = numFillRNA - numExclude;
        
        numbersData(m,:) = [i j numFillRNA numCountRNA numCyclinRNA numIntRNA];
        m = m+1;
    end
end

dlmwrite('NUMBERS_excludeNucleus.txt',numbersData,'\t');
%dlmwrite('areamatrix_excludeNucleus.txt',areaData,'\t');