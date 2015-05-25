clear dataMat;
contents = dir('data*');
tmrimgs = dir('tmr*');
cyimgs = dir('cy*');

m=1;

for i = 1:numel(contents)
    if i ==2
        continue;
    end
    load(contents(i).name);
    tmrimg = imread(tmrimgs(i).name);
    cyimg = imread(cyimgs(i).name);
    for j = 1:numel(objects)
        
        obj = objects(j);
        
        bb = obj.object_mask.boundingbox;
        
        nucmask = obj.channels.tmr.processor.mask;
        if size(nucmask,1) == 0
            continue;
        end
        
        tmrc = imcrop(tmrimg,bb);
        cyc = imcrop(cyimg,bb);
        
        tmrcrop = tmrc.*uint16(obj.object_mask.mask);
        cycrop = cyc.*uint16(obj.object_mask.mask);
        
        tmrbkgrd = min(tmrcrop(find(tmrcrop(:))));
        cybkgrd = min(cycrop(find(cycrop(:))));
        
        tmrcrop = tmrcrop-tmrbkgrd;
        cycrop = cycrop-cybkgrd;
        
        tmrnuc = tmrcrop.*uint16(nucmask);
        cynuc = cycrop.*uint16(nucmask);
        
        tmrsum = sum(tmrnuc(:));
        cysum = sum(cynuc(:));
        
        tmravg = tmrsum/numel(find(nucmask));
        cyavg = cysum/numel(find(nucmask));
        
        area = numel(find(obj.object_mask.mask==1));
        
        nucarea = numel(find(nucmask));
        
        dataMat(m,:) = [i j area nucarea tmrsum cysum tmravg cyavg];
        m = m + 1;
    end
end

dlmwrite('area_IF.txt',dataMat,'\t');