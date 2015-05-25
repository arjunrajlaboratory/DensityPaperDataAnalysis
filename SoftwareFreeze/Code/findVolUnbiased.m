function [volume height celltop_filt cellbottom_filt] = findVolUnbiased(obj,color,x,y,z);

filtersize = 30;

%x = obj.metadata.(color).x;
%y = obj.metadata.(color).y;
%z = obj.metadata.(color).z;

mask = obj.object_mask.mask;

%Find dapi mask
dapiStk = obj.channelStk('dapi');
dapiStk = max(dapiStk,[],3);
dapiMask1 = maskWithDapi(dapiStk);
dapiMask = bwareaopen(dapiMask1,2000);

points = [x y z];

R = 35;

goodpts = zeros(size(x));

for i = 1:length(x)
    currpt = points(i,:);
    for j = 1:length(x)
        dist(j) = norm(points(j,1:2)-currpt(1:2));
    end;
    ind = dist < R;  % Find all points whos XY distance is < R from currpt
    ind(i) = 0;  % Remove currpt from the list, since we aren't worried about that one.
    
    temppts = points(ind,:);
    sz = size(temppts);
    for j = 1:sz(1)
        temppts(j,:) = temppts(j,:)-currpt;
        %        temppts(j,1:2) = temppts(j,1:2)/norm(temppts(j,1:2));
    end;
    
    ind = temppts(:,3) > 0;
    if sum(ind) == 0
        goodpts(i) = 1;
    else
        
        temppts2 = temppts(ind,:);
        
        [theta,rho] = cart2pol(temppts2(:,1),temppts2(:,2));
        th = sort(theta);
        th = [th ; th(1)+2*pi];  % circularize
        maxdiff = max(diff(th));
        if maxdiff > pi
            goodpts(i) = 1;
        end;
    end;
    
    zs = temppts(:,3);  % relative Z coordinates of all neighbors (to currpt)
    if min(abs(zs)) > 5  % remove hot pixels
        goodpts(i) = 0;
    end;
    
end;

goodpts = goodpts>0;  % Convert to logical

x2 = x(goodpts);
y2 = y(goodpts);
z2 = z(goodpts);

F = TriScatteredInterp(x2,y2,z2);

ti = 1:size(mask,1);
tf = 1:size(mask,2);
[qx,qy] = meshgrid(tf,ti);
celltop = F(qx,qy);

celltop_filt = medfilt2(celltop,[filtersize filtersize]); %median filter

%mesh(qx,qy,qz);

%Now repeat for the bottom of the cell

goodpts = zeros(size(x));

for i = 1:length(x)
    currpt = points(i,:);
    for j = 1:length(x)
        dist(j) = norm(points(j,1:2)-currpt(1:2));
    end;
    ind = dist < R;  % Find all points whos XY distance is < R from currpt
    ind(i) = 0;  % Remove currpt from the list, since we aren't worried about that one.
    
    temppts = points(ind,:);
    sz = size(temppts);
    for j = 1:sz(1)
        temppts(j,:) = temppts(j,:)-currpt;
        %        temppts(j,1:2) = temppts(j,1:2)/norm(temppts(j,1:2));
    end;
    
    ind = temppts(:,3) < 0;
    if sum(ind) == 0
        goodpts(i) = 1; %keep this point if there are no lower points
    else
        
        temppts2 = temppts(ind,:);
        
        [theta,rho] = cart2pol(temppts2(:,1),temppts2(:,2));
        th = sort(theta);
        th = [th ; th(1)+2*pi];  % circularize
        maxdiff = max(diff(th));
        if maxdiff > pi
            goodpts(i) = 1;
        end;
    end;
    
    zs = temppts(:,3);  % relative Z coordinates of all neighbors (to currpt)
    if min(abs(zs)) > 5  % remove hot pixels
        goodpts(i) = 0;
    end;
    
end;

goodpts = goodpts>0;  % Convert to logical

x2 = x(goodpts);
y2 = y(goodpts);
z2 = z(goodpts);

F = TriScatteredInterp(x2,y2,z2);

ti = 1:size(mask,1);
tf = 1:size(mask,2);
[qx,qy] = meshgrid(tf,ti);
cellbottom = F(qx,qy);
cellbottom_filt = medfilt2(cellbottom,[filtersize filtersize]);

height = celltop_filt - cellbottom_filt;
height(isnan(height)) = 0;
height(~mask) = 0;
height(dapiMask) = 0;
volume = sum(height(:));