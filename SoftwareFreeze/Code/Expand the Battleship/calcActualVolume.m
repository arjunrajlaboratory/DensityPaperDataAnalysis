function [] = calcActualVolume(color)

%Calculate volume initially just to get the shape of the cell
%obj = objects(1);
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
        
        %if (i==16 & j==1)
        %    continue;
        %end
        
        fprintf('loading object %d\n',j);
        
        obj = objects(j);
        
        if isfield(obj.metadata,'volume')
            continue;
        end
        
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
        while ~(volOrig/volFake > 0.95 & volOrig/volFake < 1.05)
            if volOrig/volFake < 0.95
                break;
            end
            %Expand the volume
            expFactor = expFactor + 0.01;
            [ct ctf] = calcShell(obj,xtOrig,ytOrig,ztOrig,expFactor,30);
            [cb cbf] = calcShell(obj,xbOrig,ybOrig,zbOrig,expFactor,30);
            [xf yf zf] = fillVol(obj, color, ctf, cbf, numel(x), expFactor);
            [volFake xtf ytf ztf xbf ybf zbf ctfF cbfF] = findVol_Expand(obj, color, xf, yf, zf, expFactor);
        end
        
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
        
        volData(m,:) = [i j volume volumeWithNuc];
        m = m+1;
        
        %objects(j).metadata.volume20 = volume;
        objects(j).metadata.volume = volume;
    end
    fprintf('Saving %s\n',contents(i).name);
    save(contents(i).name,'objects');
end

warning(s);
dlmwrite('volume_expand_battleship.txt',volData,'\t');