mainDir = dir('1*');

for i = 1:numel(mainDir)
    if mainDir(i).isdir == 1
        disp(mainDir(i).name);
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
    if numel(objects) == 0
        load(dataDir(2).name);
    end
    
    if isfield(objects(1).channels,'nir')
        fillChannel = 'nir';
    else
        fillChannel = 'alexa';
    end
    
    calcVolume(fillChannel);
    
    cd('..');
end