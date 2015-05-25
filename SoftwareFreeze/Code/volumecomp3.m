%%
contents = dir('data*');

fillChannel = 'nir';
countChannel = 'alexa';

%[(datNum) (objNum) (natural) (linear) (tent)]
volumeData = zeros(50,5);

m = 1;

for i = 1:numel(contents)
    load(contents(i).name);
    for j = 1:numel(objects)
        
        if (i==7 & j==1)
            continue;
        end
        
        obj = objects(j);
        
        %Find dapi mask
        dapiStk = obj.channelStk('dapi');
        dapiStk = max(dapiStk,[],3);
        dapiMask1 = maskWithDapi(dapiStk);
        dapiMask = bwareaopen(dapiMask1,2000);
        
        x = obj.metadata.(fillChannel).x;
        y = obj.metadata.(fillChannel).y;
        z = obj.metadata.(fillChannel).z;
        
        rnacoords = [x y z];
        
        fprintf('calculating volume for object %d\n',m);
        
        [vol_natural cellheight celltop cellbottom]= tentconstruct_fromRNAspots_Olivia(dapiMask, obj.object_mask.mask, rnacoords, obj.channels.(fillChannel).maxlaplaceimage,'natural');
        [vol_linear cellheight celltop cellbottom]= tentconstruct_fromRNAspots_Olivia(dapiMask, obj.object_mask.mask, rnacoords, obj.channels.(fillChannel).maxlaplaceimage,'linear');
        [vol_tent cellheight cellvolume_delaunay cellheight_delaunay]= tentconstruct_fromRNAspots_tentinterp(obj.object_mask.mask, rnacoords);
        
        volumeData(m,:)= [i j vol_natural vol_linear vol_tent];
        
        m = m+1;

    end
end

dlmwrite('volumecomp_newMarshallCode_noBoEnf.txt',volumeData,'\t');