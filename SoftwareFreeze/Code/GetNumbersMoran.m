mainDir = dir('1*');
countChannel = 'alexa';

for i = 1:numel(mainDir)
    if mainDir(i).isdir == 1
        cd(mainDir(i).name)
    else
        cd('..');
        continue;
    end
    dataDir = dir('data*');
    if numel(dataDir) < 1
        cd('..');
        continue;
    end
    
    if ~isfield(obj.channels,'alexa')
        countChannel = 'tmr';
    else
        countChannel = 'alexa';
    end
    
    clear dat;
    n = 1;
    
    for j = 1:numel(dataDir)
        load(dataDir(j).name)
        for k = 1:numel(objects)
            obj = objects(k);
            dapiMask = obj.channels.dapi.processor.mask;
            numCountRNA = size(obj.channels.(countChannel).spotCoordinates,1);
            inNuc = numel(find(obj.channels.(countChannel).inNucleus(dapiMask)==1));
            numCytoRNA = numCountRNA - inNuc;
            dat(n,:) = [numCytoRNA inNuc];
            n = n+1;
        end
    end
    cd('..');
    dlmwrite([mainDir(i).name '.txt'],dat,'\t');
end