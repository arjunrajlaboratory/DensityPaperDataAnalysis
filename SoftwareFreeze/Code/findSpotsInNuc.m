function [numInNuc] = findSpotsInNuc(x,y,dapiMask)
%Find number of spots inside nucleus

ind = sub2ind(size(dapiMask),y,x);
numInNuc = numel(find(dapiMask(ind)==1));