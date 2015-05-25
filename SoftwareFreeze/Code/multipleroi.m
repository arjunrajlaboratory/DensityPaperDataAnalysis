%
contents = dir('data*');
fillChannel = 'nir';
countChannel = 'alexa';

m = 1;
volumeData = zeros(50,3);

for i = 1:numel(contents)
   load(contents(i).name);
   for j = 1:numel(objects)
        
       obj = objects(j);
        
       if ((i==7 | i==5) & j==1)
           continue;
       end
        
        mask = obj.object_mask.mask;
        
        %Find dapi mask
        dapiStk = obj.channelStk('dapi');
        dapiStk = max(dapiStk,[],3);
        dapiMask1 = maskWithDapi(dapiStk);
        dapiMask = bwareaopen(dapiMask1,2000);
        
        imagesc(mask);
        hold on;
        
        xp = obj.metadata.(fillChannel).x;
        yp = obj.metadata.(fillChannel).y;
        zp = obj.metadata.(fillChannel).z;
        
        scatter(xp,yp,4,zp,'filled');
        
        nRegionalMax = input('Input number of distinct local maxima:\n');
        nRoi = nRegionalMax - 1;
        
        if nRoi > 0
            cellheighttot = zeros(size(mask));
            
            regionList = zeros([size(mask) nRegionalMax]);
            
            for kk = 1:nRoi
                newROI = imfreehand;
                newMask = newROI.createMask;
                regionList(:,:,kk) = mask & newMask;
            end
            
            
            regionList(:,:,end) = mask & ~newMask;
            volList = zeros(1,nRegionalMax);
            
            %Now that we have masks for the multiple regions, find the spots
            %within each region
            
            for k = 1:size(regionList,3)
                regmask = regionList(:,:,k);
                x = round(xp);
                y = round(yp);
                z = round(zp);
                ind = sub2ind(size(regmask),y,x);
                x(find(regmask(ind)==0))=[];
                y(find(regmask(ind)==0))=[];
                z(find(regmask(ind)==0))=[];
                %figure; imagesc(regmask); hold on;
                %scatter(x,y,4,z,'filled');
                rnacoords = [x y z];
                [vol_natural cellheight celltop cellbottom]= tentconstruct_fromRNAspots_Olivia(dapiMask, regmask, rnacoords, obj.channels.(fillChannel).maxlaplaceimage,'natural');
                volList(k) = vol_natural;
                cellheighttot = cellheighttot + cellheight;
            end
            vol_natural = sum(volList);
        else
            rnacoords = [xp yp zp];
            [vol_natural cellheight celltop cellbottom]= tentconstruct_fromRNAspots_Olivia(dapiMask, obj.object_mask.mask, rnacoords, obj.channels.(fillChannel).maxlaplaceimage,'natural');
        end
        
        close;
        
       volumeData(m,:) = [i j vol_natural];
       m = m+1;
   end
end

dlmwrite('volume_noboenf_multroi.txt',volumeData,'\t');