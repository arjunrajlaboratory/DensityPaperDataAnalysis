function [zDiff] = findZdifference(data);
[i,j,k] = data.getSpotCoordinates();
zDiff = max(k)-min(k);
end