contents = dir('data*');
fillChannel = 'nir';
countChannel = 'alexa';

volumeData = zeros(50,3);
n = 1;

for k = 1:numel(contents)
    load(contents(k).name);
    for kk = 1:numel(objects)
        
        fprintf('calculating volume for object %d\n',n);
        
        obj = objects(kk);
        
        if ~isfield(obj.metadata,fillChannel)
            continue;
        end
        
        mask = obj.object_mask.mask;
        
        %Find dapi mask
        dapiStk = obj.channelStk('dapi');
        dapiStk = max(dapiStk,[],3);
        dapiMask1 = maskWithDapi(dapiStk);
        dapiMask = bwareaopen(dapiMask1,2000);
        
        x = obj.metadata.(fillChannel).x;
        y = obj.metadata.(fillChannel).y;
        z = obj.metadata.(fillChannel).z;
        
        xdim = size(mask,2);
        ydim = size(mask,1);
        
        step = 20;
        
        numBoxesInMask = 0;
        vol = 0;
        height = zeros(size(mask));
        %mins = 0;
        %n = 0;
        
        for i = 1:step:xdim
            for j = 1:step:ydim
                minimask = zeros(size(mask));
                imax = min(i+step-1,xdim);
                jmax = min(j+step-1,ydim);
                %minimask = poly2mask(i,j,imax,jmax);
                minimask(j:jmax,i:imax) = 1;
                if ~(mask & minimask)
                    continue;
                end
                numBoxesInMask = numBoxesInMask + 1;
                inpts = inpolygon(x,y,[i imax],[j jmax]);
                
                %**this is an extra smoothing step**
                if isempty(find(inpts)) | max(z(inpts)) - min(z(inpts)) == 0
                    %expand box
                    imin = max(1,i-step);
                    imax = min(imax+step,xdim);
                    jmin = max(1,j-step);
                    jmax = min(jmax+step,ydim);
                    inpts = inpolygon(x,y,[imin imax],[jmin jmax]);
                end
                if isempty(find(inpts))
                    continue;
                else
                    %mins = mins + min(z(inpts));
                    %n = n + 1;
                    height(mask & minimask) = max(z(inpts)) - min(z(inpts));
                end
            end
        end
        height(dapiMask) = 0;
        volumeData(n,:) = [k kk sum(height(:))];
        n = n+1;
    end
end

dlmwrite('volume_grid.txt',volumeData,'\t');