mainDir = dir('1*');

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
    
    load(dataDir(1).name);
    if numel(objects) < 1
        load(dataDir(2).name);
    end
    
    if isfield(objects(1).channels,'alexa')
        countChannel = 'alexa';
        intChannel = 'tmr';
    else
        countChannel = 'nir';
        intChannel = 'tmr';
    end
    
    % Check to see that data has been processed and reviewed
    if objects(1).channels.(countChannel).isProcessed & objects(1).channels.(countChannel).processor.reviewed & ~isfield(objects(1).channels.(countChannel).metadata,'gaussFitPostProc')
        postprocessimageobjects('channels',{countChannel,intChannel},'imageProcessors',{'gaussFitPostProc','gaussFitPostProc'},'GUI',false);
    end
    cd('..');
end