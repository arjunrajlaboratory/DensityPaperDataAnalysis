function [] = BinarySearchVol(color);

%color = 'nir';
contents = dir('data*');
m = 1;
clear volData;

msgid = 'MATLAB:TriScatteredInterp:DupPtsAvValuesWarnId';
s = warning('off',msgid);

for i = 1:numel(contents)
    load(contents(i).name);
    fprintf('loading %s\n',contents(i).name);
    for j = 1:numel(objects)
        flag = 0;
        
        %if (i == 19 | i == 20 | i == 24) & j == 1
        %    continue;
        %end
        
        fprintf('loading object %d\n',j);
        
        obj = objects(j);
        
        if isfield(obj.metadata,'volumebinary')
            continue;
        end
        
        ts = tic;
        
        x = obj.channels.(color).spotCoordinates(:,2);
        y = obj.channels.(color).spotCoordinates(:,1);
        z = obj.channels.(color).spotCoordinates(:,3);
        
        [volOrig xtOrig ytOrig ztOrig xbOrig ybOrig zbOrig ctfOrig cbfOrig] = findVol_Expand(obj, color, x, y, z, 1);
        
        %Randomly fill the volume with different number of spots until the number
        %of "fake" spots is within 10% of the number of actual spots
        [xf yf zf] = fillVol(obj, color, ctfOrig, cbfOrig, numel(x), 1);
        
        %We want "real" (original) volume / "fake" volume ratio -> 1 (say to within 5%)
        
        %Now that we have the right xf, yf, zf, calculate volume
        expFactor = 1;
        [volFake xtf ytf ztf xbf ybf zbf ctf cbf] = findVol_Expand(obj, color, xf, yf, zf, expFactor);
        
        %Expand the cell until the two volumes are close to each other
        
        %find the end of the interval, start with expFactor = 2*begInt
        begInt = expFactor;
        endInt = 1.1*begInt;
        while 'true'
            [volFake xtf ytf ztf xbf ybf zbf ctf cbf] = findVol_Expand(obj, color, xf, yf, zf, endInt);
            if abs(1-volOrig/volFake)<0.01 % this means we've found the right volume
                flag = 1;
                break;
            elseif volOrig/volFake < 1 %this means endInt is too large, so keep as end
                break;
            else %endInt is too small, keep expanding
                endInt = endInt+0.1;
            end
        end
        
        %now we have an endInt and a begInt
        tic;
        while 'true'
            if flag == 1
                break
            end
            expFactor = (endInt+begInt)/2;
            [volFake xtf ytf ztf xbf ybf zbf ctf cbf] = findVol_Expand(obj, color, xf, yf, zf, expFactor);
            if abs(1-volOrig/volFake)<0.01
                break;
            elseif volOrig/volFake > 1 %this means our expFactor is too small
                begInt = expFactor; %make this the beginning of our interval
                expFactor = (endInt+begInt)/2; %increase expFactor
            else %this means our expFactor is too large
                endInt = expFactor; %make this the end of the interval
                expFactor = (endInt+begInt)/2;
            end
        end
        toc;
        
        %Now the actual volume is what is defined by the most current celltop and
        %cellbottom maps
        
        mask = imresize(obj.object_mask.mask,expFactor);
        
        %Find dapi mask
%         dapiStk = obj.channelStk('dapi');
%         dapiStk = max(dapiStk,[],3);
%         dapiMask1 = maskWithDapi(dapiStk);
%         dapiMask = bwareaopen(dapiMask1,2000);
%         dapiMask = imresize(dapiMask,expFactor);
        dapiMask = imresize(obj.channels.dapi.processor.mask,expFactor);
        
        height = ctf - cbf;
        height(isnan(height)) = 0;
        height(~mask) = 0;
        volumeWithNuc = sum(height(:));
        height(dapiMask) = 0;
        volume = sum(height(:));
        
        toc(ts);
        
        volData(m,:) = [i j volume volumeWithNuc];
        m = m+1;
        
        %objects(j).metadata.volume20 = volume;
        objects(j).metadata.volumebinary = volume;
    end
    fprintf('Saving %s\n',contents(i).name);
    save(contents(i).name,'objects');
end

warning(s);
dlmwrite('volume_expand_battleship.txt',volData,'\t');
