function [xf yf zf] = genUnifPtsExact(obj, color, ctf, cbf, nSpots, expFactor)

% Expand mask
mask = imresize(obj.object_mask.mask,expFactor);

%Generate random points filling volume xdim * ydim * max(z)

%this corresponds to maximum x value
xdim = size(mask,2)-1;

%this corresponds to max y value
ydim = size(mask,1)-1;

zdim = round(max(obj.channels.(color).spotCoordinates(:,3))*expFactor)-1;
n = 1;

while n < nSpots
    
    %Add a point
    xn = xdim*rand()+1;
    yn = ydim*rand()+1;
    
    %Now mask by x,y
    xr = round(xn);
    yr = round(yn);
    
    %Check if the point is in the mask
    ind = sub2ind(size(mask),yr,xr);
    if mask(ind)~=1
        continue;
    end
    
    %Now we need to select only z points which are of the correct height
    zn = zdim*rand()+1;
    
    if zn <= ctf(yr,xr) && zn >= cbf(yr,xr)
        xf(n,:) = xn;
        yf(n,:) = yn;
        zf(n,:) = zn;
        n = n + 1;
    end
end