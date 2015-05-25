function [ cellvolume cellheight cellvolume_delaunay cellheight_delaunay] = tentconstruct_fromRNAspots_tentinterp( objmask , rnacoords )
% obj is the image object representing for example a cell. 
% rnacoords is an Mx3 with the [x y z] coordinates of every RNA spot to be
% used in constructing the top and bottom envelopes (tents) around the
% object. 

% The difference between this function and the old
% tentconstruct_fromRNAspots is that, in addition to the Delaunay
% triangulation, it does a traingulation that is true to the tent
% construction algorithm. The cell volume and height calculated with the
% new triangulation are output as cellvolume and cellheight.  The Delaunay
% based interpolation results are output as cellvolume_delaunay and
% cellheight_delaunay, but we use the 'linear', 
% not the 'natural' interpolation method. 

xp=rnacoords(:,1)';
yp=rnacoords(:,2)';
zp=rnacoords(:,3)';

% First create the tent at the top.

% The next line is all that needs to change to do the bottom tent
% calculation. 'descend'-top tent, 'ascend' - bottom tent. 
[zp_sorted,ind_sort]=sort(zp,'descend');
xp_sorted=xp(ind_sort);
yp_sorted=yp(ind_sort);

% 'out' subscript means points that could possibly become a vertex of the
% tent.
zp_out=zp_sorted(4:end);
xp_out=xp_sorted(4:end);
yp_out=yp_sorted(4:end);

% 'tent' subscript are the points that make up the tent. Start with the
% three highest points.
zp_tent=zp_sorted(1:3);
xp_tent=xp_sorted(1:3);
yp_tent=yp_sorted(1:3);

% the first triangle contains just the first three points in the tent.
tenttriangulation=[1 2 3];
hullindicesold=convhull(xp_tent,yp_tent);

while ~isempty(zp_out)
    
    % find all points not in the tent that are within its convex hull.
    hullindices=convhull(xp_tent,yp_tent);
    out_within_tent=inpolygon(xp_out,yp_out,...
        xp_tent(hullindices),yp_tent(hullindices));
    % eliminate those points from future consideration
    zp_out(out_within_tent)=[];
    xp_out(out_within_tent)=[];
    yp_out(out_within_tent)=[];
    
    % now determine the new triangles that have been created. do not do
    % this step on the first iteration of this loop (when there are three
    % in the tent). 
    if length(zp_tent)>3
        newest_tentpoint=length(zp_tent);
        tenttriangulation=[tenttriangulation; ...
            findnewtriangles_tent(hullindices,hullindicesold,newest_tentpoint)];
    end
    
    if ~isempty(zp_out); 
        % keep the old convex hull
        hullindicesold=hullindices;
        % Add the highest point in the 'out' set. to the tent.
        zp_tent=[zp_tent zp_out(1)];
        xp_tent=[xp_tent xp_out(1)];
        yp_tent=[yp_tent yp_out(1)];
        % remove it from the 'out' set. 
        zp_out(1)=[];
        xp_out(1)=[];
        yp_out(1)=[];
    end
end
%Store the top tent and the Delaunay triangulation of it for future
%display. 
toptent=[xp_tent' yp_tent' zp_tent'];
topdt=DelaunayTri(xp_tent',yp_tent');
toptenttriangle=tenttriangulation;



% ------------   Repeating for BOTTOM tent -----

% The next line is all that needs to change to do the bottom tent
% calculation. 'descend'-top tent, 'ascend' - bottom tent. 
[zp_sorted,ind_sort]=sort(zp,'ascend');
xp_sorted=xp(ind_sort);
yp_sorted=yp(ind_sort);

% 'out' subscript means points that could possibly become a vertex of the
% tent.
zp_out=zp_sorted(4:end);
xp_out=xp_sorted(4:end);
yp_out=yp_sorted(4:end);

% 'tent' subscript are the points that make up the tent. Start with the
% three highest points.
zp_tent=zp_sorted(1:3);
xp_tent=xp_sorted(1:3);
yp_tent=yp_sorted(1:3);

% the first triangle contains just the first three points in the tent.
tenttriangulation=[1 2 3];
hullindicesold=convhull(xp_tent,yp_tent);

while ~isempty(zp_out)
    
    % find all points not in the tent that are within its convex hull.
    hullindices=convhull(xp_tent,yp_tent);
    out_within_tent=inpolygon(xp_out,yp_out,...
        xp_tent(hullindices),yp_tent(hullindices));
    % eliminate those points from future consideration
    zp_out(out_within_tent)=[];
    xp_out(out_within_tent)=[];
    yp_out(out_within_tent)=[];
    
    % now determine the new triangles that have been created. do not do
    % this step on the first iteration of this loop (when there are three
    % in the tent). 
    if length(zp_tent)>3
        newest_tentpoint=length(zp_tent);
        tenttriangulation=[tenttriangulation; ...
            findnewtriangles_tent(hullindices,hullindicesold,newest_tentpoint)];
    end
    
    if ~isempty(zp_out); 
        % keep the old convex hull
        hullindicesold=hullindices;
        % Add the highest point in the 'out' set. to the tent.
        zp_tent=[zp_tent zp_out(1)];
        xp_tent=[xp_tent xp_out(1)];
        yp_tent=[yp_tent yp_out(1)];
        % remove it from the 'out' set. 
        zp_out(1)=[];
        xp_out(1)=[];
        yp_out(1)=[];
    end
end

%Store the top tent and the Delaunay triangulation of it for future
%display. 
bottomtent=[xp_tent' yp_tent' zp_tent'];
bottomdt=DelaunayTri(xp_tent',yp_tent');
bottomtenttriangle=tenttriangulation;




% the volume is calculated by building interpolation functions for
% both the top and the bottom and creating pixel by pixel maps of the top 
% and bottom tents.

% Using the 'natural' method because it gives a much smoother
% interpolation. This is the Delaunay triangulation based method. We'll do
% the tent-based triangulation and interpolation after. 
% Ftop=TriScatteredInterp(topdt,toptent(:,3),'natural');
% Fbottom=TriScatteredInterp(bottomdt,bottomtent(:,3),'natural');
Ftop=TriScatteredInterp(topdt,toptent(:,3),'linear');
Fbottom=TriScatteredInterp(bottomdt,bottomtent(:,3),'linear');

Xs=repmat(1:size(objmask,2),[size(objmask,1) 1]);
Ys=repmat((1:size(objmask,1))',[1 size(objmask,2)]);

% interpolate the tent heights at every pixel of the cell.
celltop_delaunay=Ftop(Xs,Ys);
cellbottom_delaunay=Fbottom(Xs,Ys);
% the interpolation procedure produces NaN for points outside the
% triangulation, so we set those to zero instead:
celltop_delaunay(isnan(celltop_delaunay))=0;
cellbottom_delaunay(isnan(cellbottom_delaunay))=0;

cellheight_delaunay=celltop_delaunay-cellbottom_delaunay;
cellheight_delaunay(~(objmask))=0;

cellvolume_delaunay=sum(cellheight_delaunay(:));


% Now do the tent based triangulation.
% Start with a NaN matrix
celltop=nan(size(Xs));

% Loop through all the triangles in the triangulation
for trianglenum=1:size(toptenttriangle,1)
    % linear indices of all points with no value assigned. 
    unassignedpoints=find(isnan(celltop(:)));
    
    % Coordinates of current triangle
    trianglexcoords=toptent(toptenttriangle(trianglenum,:),1);
    triangleycoords=toptent(toptenttriangle(trianglenum,:),2);
    trianglezcoords=toptent(toptenttriangle(trianglenum,:),3);
    % Interpolation function for this triangle
    Fsingletriangle=TriScatteredInterp(trianglexcoords,triangleycoords,...
        trianglezcoords);
    % Use triangle to assign some points that have not been assigned yet.
    celltop(unassignedpoints)=Fsingletriangle(Xs(unassignedpoints),...
        Ys(unassignedpoints));
end

% Now do the tent based triangulation.
% Start with a NaN matrix
cellbottom=nan(size(Xs));

% Loop through all the triangles in the triangulation
for trianglenum=1:size(bottomtenttriangle,1)
    % linear indices of all points with no value assigned. 
    unassignedpoints=find(isnan(cellbottom(:)));
    
    % Coordinates of current triangle
    trianglexcoords=bottomtent(bottomtenttriangle(trianglenum,:),1);
    triangleycoords=bottomtent(bottomtenttriangle(trianglenum,:),2);
    trianglezcoords=bottomtent(bottomtenttriangle(trianglenum,:),3);
    % Interpolation function for this triangle
    Fsingletriangle=TriScatteredInterp(trianglexcoords,triangleycoords,...
        trianglezcoords);
    % Use triangle to assign some points that have not been assigned yet.
    cellbottom(unassignedpoints)=Fsingletriangle(Xs(unassignedpoints),...
        Ys(unassignedpoints));
end

celltop(isnan(celltop))=0;
cellbottom(isnan(cellbottom))=0;

cellheight=celltop-cellbottom;
cellheight(~(objmask))=0;

% figure(4);
% subplot(1,3,1);imshow(celltop/max(celltop(:)));
% subplot(1,3,2);imshow(cellbottom/max(celltop(:)));
% subplot(1,3,3);imshow(cellheight/max(celltop(:)));
% colormap jet;

cellvolume=sum(cellheight(:));

end

