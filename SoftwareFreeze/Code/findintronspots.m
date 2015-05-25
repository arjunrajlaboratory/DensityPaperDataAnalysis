%load data014;
%obj = objects(1);

dapiMask = obj.channels.dapi.processor.mask;

countChannel = 'alexa';
intChannel = 'tmr';

%Threshold for whether spots are colocalized
thresh = 50;


%All mRNA spots
mCoords = obj.channels.(countChannel).spotCoordinates;
%Choose only those in nucleus
mrc = round(mCoords);
x = mrc(:,2);
y = mrc(:,1);
ind = sub2ind(size(dapiMask),y,x);
mCoords(find(dapiMask(ind)==0),:) = [];

%All intron spots
iCoords = obj.channels.(intChannel).spotCoordinates;

%Scale
icp = iCoords;
icp(:,3) = icp(:,3)*.35/.125;

%Choose only those in nucleus
irc = round(icp);
x = irc(:,2);
y = irc(:,1);
ind = sub2ind(size(dapiMask),y,x);
icp(find(dapiMask(ind)==0),:) = [];
iCoords(find(dapiMask(ind)==0),:) = [];

%Binary list of "real" intron spots
keep = false(size(iCoords,1),1);

%Find nearest neighbor mRNA spot for each intron spot
nnIdx = knnsearch(mCoords,icp);

%Keep the spots that are colocalized
for i = 1:size(iCoords,1)
    ic = icp(i,:);
    mc = mCoords(nnIdx(i),:);
    if pdist([ic ; mc]) < thresh
        keep(i) = true;
    end
end

iCoords = iCoords(keep,:);

%%

contents = dir('data*');
m = 1;
clear numIntList

for i = 1:numel(contents)
    load(contents(i).name);
    for j = 1:numel(objects)
        
        obj = objects(j);
        
        %Returns indices of tmr spots in nucleus
        dapiMask = obj.channels.dapi.processor.mask;
        a1 = obj.channels.tmr.inNucleus(dapiMask);
        
        %Returns indices of tmr spots colocalized with alexa spots
        a2 = obj.channels.tmr.getColocalized(obj.channels.alexa,3,.35/.125);
        
        %Intron spot coordinates
        iCoords = obj.channels.tmr.spotCoordinates(a1&a2,:);
        
        numIntList(m,:) = [i j size(iCoords,1)];
        m = m+1;
    end
end

dlmwrite('numInt.txt',numIntList,'\t');
        
        
