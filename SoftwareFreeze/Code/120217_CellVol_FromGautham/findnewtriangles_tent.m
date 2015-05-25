function [newtriangles] = findnewtriangles_tent(newhull_raw,oldhull_raw,newpoint)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% suppose we started with N points in (x,y), and found their convexhull
% using the convhull function, such that  oldhull_raw=convhull(X,Y), and we
% then add another point that is outside of the oldhull, so that there are 
% N+1 points, and there is a new convexhull newhull_raw=convhull([X xnew],[Y ynew]).
% We seek to find the vertex indices of all the new triangles that have to be
% built to fill in the newhull starting from the old hull. 
% The new point should be = N+1.

% the algorithm exploits the fact that convhull returns the indices of the
% vertices making up the edge of the hull in counterclockwise order. 

% Suppose that newhull has the structure:
%  newhull =  ... a , newpoint, b, ... 
%  and that the oldhull has the structure:
%  oldhull = ... a , c_1, ..., c_k, b, ...
%  then we need to add the triangles:
%   [a, c_1, newpoint]
%   [c_1, c_2, newpoint]
%   ...
%   [c_(k-1), c_k, newpoint]
%   [c_k, b, newpoint]
%  which will be a (k+1, 3) size matrix. 

newtriangles=zeros(0,3);

% Remove the duplicate entry in newhull and oldhull: (the first entry is
% always duplicated at the end). 
newhull=newhull_raw(1:end-1);
oldhull=oldhull_raw(1:end-1);

% First find the newpoint in the newhull:

newpointpos_innewhull=find(newhull==newpoint);

if isempty(newpointpos_innewhull)
    return
end

%find the values of a and b: Use modular arythmetic to handle cyclical
%behavior of newhull, oldhull;
prevpointpos_innewhull=1+mod(newpointpos_innewhull-1-1,length(newhull));
prevpoint=newhull(prevpointpos_innewhull);
nextpointpos_innewhull=1+mod(newpointpos_innewhull+1-1,length(newhull));
nextpoint=newhull(nextpointpos_innewhull);

% Now find the positions of a and b in the old hull:

prevpointpos_inoldhull=find(oldhull==prevpoint);
nextpointpos_inoldhull=find(oldhull==nextpoint);

if isempty(prevpointpos_inoldhull)||isempty(nextpointpos_inoldhull)
    return
end 

% in case b is looped around and finds itself before a in the old hull
if nextpointpos_inoldhull<prevpointpos_inoldhull
   nextpointpos_inoldhull=nextpointpos_inoldhull+length(oldhull); 
end

% The positions in the old hull of [a, c_1, ..., c_k, b]
oldhullpos_totriangulate=prevpointpos_inoldhull:nextpointpos_inoldhull;
oldhullpos_totriangulate=1+mod(oldhullpos_totriangulate-1,length(oldhull));


for i=1:length(oldhullpos_totriangulate)-1
    newtriangles = [newtriangles; ...
        oldhull(oldhullpos_totriangulate(i)) ...
        oldhull(oldhullpos_totriangulate(i+1)) newpoint];
end





end

