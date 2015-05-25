mainDir = dir('1*');
clear dat;
dat = {};

for i = 1:numel(mainDir)
    if mainDir(i).isdir == 1
        cd(mainDir(i).name)
        disp(mainDir(i).name)
    else
        cd('..');
        continue;
    end
    dataDir = dir('data*');
    if numel(dataDir) < 1
        cd('..');
        continue;
    end
    
    st = regexp(mainDir(i).name,'_','split');
    
    load(dataDir(1).name);
    if numel(fieldnames(objects(1).channels)) == 7
        countChannel = 'alexa';
        cyclinChannel = 'cy';
        intChannel = 'tmr';
    else
        if strcmp(st(2),'GAPDH')
            countChannel = 'nir';
            cyclinChannel = 'cy';
            intChannel = 'tmr';
        elseif strcmp(st(2),'POLR2A')
            cd('..');
            continue;
        end
    end
    
    [fp,efp,fn,efn,ip,eip,in,ein]=Frequency_Intensity_AnyGene(countChannel,cyclinChannel,intChannel);
    dat{i,1} = st(2);
    dat{i,2} = fp;
    dat{i,3} = efp;
    dat{i,4} = fn;
    dat{i,5} = efn;
    dat{i,6} = ip;
    dat{i,7} = eip;
    dat{i,8} = in;
    dat{i,9} = ein;
    
    cd('..');
    %dlmwrite([mainDir(i).name '.txt'],dat,'\t');
end