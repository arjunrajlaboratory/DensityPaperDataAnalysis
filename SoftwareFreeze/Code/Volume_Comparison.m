%Volume comparison
% -Natural Delaunay
% -Natural Delaunay with boundary enforced
% -Linear Delaunay
% -Linear Delaunay with boundary enforced
% -Tent-based (here, cannot enforce boundary)

contents = dir('data*');
rna_color='cy';

%Find the total number of objects in this folder
n = 0;
for i = 1:numel(contents)
    load(contents(i).name);
    n = n + numel(objects);
end

volumeComp = zeros(n,4);

m = 1;

%Make a matrix with volume/number/density data
%for i = 1:numel(contents)
for i = 1:5
    load(contents(i).name);
    for j = 1:numel(objects)
        
        obj=objects(j);
        
        dapiStk = obj.channelStk('dapi');
        dapiStk = max(dapiStk,[],3);
        dapiMask1 = maskWithDapi(dapiStk);
        dapiMask = bwareaopen(dapiMask1,2000);
        
        
        xp=obj.channels.(rna_color).fitdataRNAonly.yp_fit;
        yp=obj.channels.(rna_color).fitdataRNAonly.xp_fit;
        zp=obj.channels.(rna_color).fitdataRNAonly.rawzp;
        
        rnacoords=[xp', yp', zp'];
        %[cellvolume_tent cellheight_tent cellvolume_delaunaylinear cellheight_delaunaylinear]= tentconstruct_fromRNAspots_tentinterp( obj.object_mask.mask , rnacoords );
        %[cellvolume_delaunaynatural cellheight_delaunaynatural]= tentconstruct_fromRNAspots( obj.object_mask.mask , rnacoords );
        [cellvolume_boundarynatural cellheight_boundarynatural] = objectvolume_fromRNAspots_Olivia(obj,rna_color,dapiMask);
        cellvolume_boundarynatural
        %volumeComp(m,:)=[cellvolume_delaunaynatural cellvolume_boundarynatural cellvolume_delaunaylinear cellvolume_tent];
        %fprintf('Calculating density for object %d of %d\n',m,n);
        m = m+1;
    end
end