clear;
load data014

obj = objects(1);
fillChannel = 'nir';
countChannel = 'tmr';

n = 20;

xp = obj.metadata.(fillChannel).x;
yp = obj.metadata.(fillChannel).y;
zp = obj.metadata.(fillChannel).z;


mask = obj.object_mask.mask;

%Find dapi mask
dapiStk = obj.channelStk('dapi');
dapiStk = max(dapiStk,[],3);
dapiMask1 = maskWithDapi(dapiStk);
dapiMask = bwareaopen(dapiMask1,2000);

mask(dapiMask) = 0;

for i = 1:n
    
    imagesc(mask);
    newROI = imfreehand;
    newMask = newROI.createMask;
    maskList(:,:,i) = mask & newMask;
    close;
end

for i = 1:n
    thisMask = maskList(:,:,i);
    x = round(xp);
    y = round(yp);
    z = round(zp);
    ind = sub2ind(size(thisMask),y,x);
    x(find(thisMask(ind)==0))=[];
    y(find(thisMask(ind)==0))=[];
    z(find(thisMask(ind)==0))=[];
    %figure; imagesc(regmask); hold on;
    %scatter(x,y,4,z,'filled');
    numRNA = numel(x);
    [volume height] = findVolUnbiased(obj,fillChannel);
    height(~thisMask) = 0;
    volume = sum(height(:));
    decVol(i,:) = [volume numRNA]
end

%plot

[b, bint] = regress(decVol(:,2), [ones(length(decVol),1), decVol(:,1)], 0.05);
scatter(decVol(:,1),decVol(:,2));
hold on;
refline(b(2),b(1))
h1 = refline(bint(2,1),bint(1,2));
h2 = refline(bint(2,2),bint(1,1));
set(h1,'Color','r')
set(h2,'Color','r')