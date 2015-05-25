%% 
   %This counts the number of RNA using Marshall's new method, then
   %calculates volume using a particular method.
   %The object's metadata is changed to reflect the positions of the RNA
   %that has been counted.

contents = dir('data*');

fillChannel = 'nir';
countChannel = 'alexa';

%Find the total number of objects in this folder
n = 0;
for i = 1:numel(contents)
    load(contents(i).name);
    n = n + numel(objects);
end

%volumeData = zeros(n,7);

%[(datNum) (objNum) (volume) (fillCount) (countCount)]
volumeData = zeros(n,5);

m = 1;

for i = 1:numel(contents)
    load(contents(i).name);
    for j = 1:numel(objects)
        
        if (i==7 & j==1)
            continue;
        end
        
        obj = objects(j);
        
        %Find spot locations of 'fill' channel
        [nspots x y z] = findSpotsATrousCC(obj,fillChannel,'lenient');
        obj.metadata.(fillChannel).x = x;
        obj.metadata.(fillChannel).y = y;
        obj.metadata.(fillChannel).z = z;
        
        %Find dapi mask
        dapiStk = obj.channelStk('dapi');
        dapiStk = max(dapiStk,[],3);
        dapiMask1 = maskWithDapi(dapiStk);
        dapiMask = bwareaopen(dapiMask1,2000);
        
        rnacoords = [x y z];
        
        [cellvolume cellheight celltop cellbottom]= tentconstruct_fromRNAspots_Olivia(dapiMask, obj.object_mask.mask, rnacoords, obj.channels.(fillChannel).maxlaplaceimage);
        %[cellvolume cellheight celltop cellbottom]= tentconstruct_fromRNAspots_forceOutline_Olivia(dapiMask, obj.object_mask.mask, rnacoords, obj.channels.(fillChannel).maxlaplaceimage);
        
        [nspots x y z] = findSpotsATrousCC(obj,countChannel,'lenient');
        obj.metadata.(countChannel).x = x;
        obj.metadata.(countChannel).y = y;
        obj.metadata.(countChannel).z = z;
        
        numFillRNA = numel(obj.metadata.(fillChannel).x);
        numCountRNA = numel(obj.metadata.(countChannel).x);
        
        x = round(obj.metadata.(countChannel).y);
        y = round(obj.metadata.(countChannel).x);
        
        xD = round(obj.metadata.(fillChannel).y);
        yD = round(obj.metadata.(fillChannel).x);
        
        
        ind = sub2ind(size(dapiMask),x,y);
        numExclude = numel(find(dapiMask(ind)==1));
        
        numCountRNA = numCountRNA - numExclude;
        
        ind = sub2ind(size(dapiMask),xD,yD);
        numExclude = numel(find(dapiMask(ind)==1));
        
        numFillRNA = numFillRNA - numExclude;
        
        volumeData(m,:)= [i j cellvolume numFillRNA numCountRNA];
        fprintf('Calculating density for object %d of %d\n',m,n);
        m = m+1;
    end
end

dlmwrite('densitymatrix_newMarshallCode_noBoEnf.txt',volumeData,'\t');