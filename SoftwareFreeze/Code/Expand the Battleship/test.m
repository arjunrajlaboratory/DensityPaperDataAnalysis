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
tic;
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
toc;

%Now the actual volume is what is defined by the most current celltop and
%cellbottom maps

mask = imresize(obj.object_mask.mask,expFactor);

%Find dapi mask
dapiStk = obj.channelStk('dapi');
dapiStk = max(dapiStk,[],3);
dapiMask1 = maskWithDapi(dapiStk);
dapiMask = bwareaopen(dapiMask1,2000);
dapiMask = imresize(dapiMask,expFactor);

height = ctf - cbf;
height(isnan(height)) = 0;
height(~mask) = 0;
volumeWithNuc = sum(height(:));
height(dapiMask) = 0;
volume = sum(height(:));

toc(ts);