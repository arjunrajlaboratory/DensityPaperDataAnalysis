function [volume x2t y2t z2t x2b y2b z2b celltop cellbottom goodpts_top goodpts_bottom] = findVol_Expand(obj, color, x, y, z, expFactor);

%Returns (volume), (coords on outside of "top" shell), (coords on outside
%of "bottom" shell), (heatmap of top of cell), (heatmap of bottom)

filtersize = 30;

%x = obj.metadata.(color).x;
%y = obj.metadata.(color).y;
%z = obj.metadata.(color).z;

mask = imresize(obj.object_mask.mask,expFactor);

%Find dapi mask
dapiStk = obj.channelStk('dapi');
dapiStk = max(dapiStk,[],3);
dapiMask1 = maskWithDapi(dapiStk);
dapiMask = bwareaopen(dapiMask1,2000);
dapiMask = imresize(dapiMask,expFactor);

points = [x y z];

%R = 35 %changed to 25 8/20/12
R = 20;

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

x2t = x(goodpts);
y2t = y(goodpts);
z2t = z(goodpts);

goodpts_top = goodpts;

% x2 = x2*expFactor;
% y2 = y2*expFactor;
% z2 = z2*expFactor;
% 
% F = TriScatteredInterp(x2,y2,z2);
% 
% ti = 1:ceil(size(mask,1)*expFactor);
% tf = 1:ceil(size(mask,2)*expFactor);
% [qx,qy] = meshgrid(tf,ti);
% celltop = F(qx,qy);
% celltop_filt = medfilt2(celltop,[filtersize filtersize]); %median filter

[celltop celltop_filt] = calcShell(obj,x2t,y2t,z2t,expFactor,filtersize);

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

x2b = x(goodpts);
y2b = y(goodpts);
z2b = z(goodpts);

goodpts_bottom = goodpts;

% x2 = x2*expFactor;
% y2 = y2*expFactor;
% z2 = z2*expFactor;
% 
% F = TriScatteredInterp(x2,y2,z2);
% 
% ti = 1:ceil(size(mask,1)*expFactor);
% tf = 1:ceil(size(mask,2)*expFactor);
% [qx,qy] = meshgrid(tf,ti);
% cellbottom = F(qx,qy);
% cellbottom_filt = medfilt2(cellbottom,[filtersize filtersize]);

[cellbottom cellbottom_filt] = calcShell(obj,x2b,y2b,z2b,expFactor,filtersize);

height = celltop_filt - cellbottom_filt;
height(isnan(height)) = 0;
height(~mask) = 0;
height(dapiMask) = 0;
volume = sum(height(:));

%if nargout > 1
%    varargout = [volume x2t y2t z2t x2b y2b z2b celltop_filt cellbottom_filt];
%end