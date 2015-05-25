%This is the objectvolume_fromRNAspots/tentconstruct_fromRNAspots script
%without functions.

obj = objects(1);
rna_color = 'nir';

xp=obj.channels.(rna_color).fitdataRNAonly.yp_fit;
yp=obj.channels.(rna_color).fitdataRNAonly.xp_fit;
zp=obj.channels.(rna_color).fitdataRNAonly.rawzp;

rnacoords=[xp', yp', zp'];

objmask = obj.object_mask.mask;

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

while ~isempty(zp_out)
    
    % find all points not in the tent that are within the tent's convex hull
    
    hullindices=convhull(xp_tent,yp_tent);
    out_within_tent=inpolygon(xp_out,yp_out,...
        xp_tent(hullindices),yp_tent(hullindices));
    
    % eliminate those points from future consideration
    
    zp_out(out_within_tent)=[];
    xp_out(out_within_tent)=[];
    yp_out(out_within_tent)=[];
    
    
    if ~isempty(zp_out);
        
        % Add the highest point in the 'out' set. to the tent.
        zp_tent=[zp_tent zp_out(1)];
        xp_tent=[xp_tent xp_out(1)];
        yp_tent=[yp_tent yp_out(1)];
        
        % remove it from the 'out' set. 
        zp_out(1)=[];
        xp_out(1)=[];
        yp_out(1)=[];
    end
    
    % keep going until all points are in tent or excluded
end

%Store the top tent and the Delaunay triangulation of it for future
%display. 

%Find mask outline and force points on the outline to be included
%in the Delaunay calculation

B = bwboundaries(objmask);
b1 = B{1};
x = b1(:,2);
y = b1(:,1);
x1 = x(1:20:end);
y1 = y(1:20:end);
r = repmat(min(zp_tent), numel(x1), 1);

xp_tent = [xp_tent x1'];
yp_tent = [yp_tent y1'];
zp_tent = [zp_tent r'];

toptent=[xp_tent' yp_tent' zp_tent'];
topdt=DelaunayTri(xp_tent',yp_tent');



% ------------   Repeating for BOTTOM tent -----

%Switching this line to ascend and repeating everything gives the bottom
%tent.
[zp_sorted,ind_sort]=sort(zp,'ascend');
xp_sorted=xp(ind_sort);
yp_sorted=yp(ind_sort);

% 'out' subscript means points that could possibly become a vertex of the
% tent.
zp_out=zp_sorted(4:end);
xp_out=xp_sorted(4:end);
yp_out=yp_sorted(4:end);

% 'tent' subscript are the points that make up the tent. Start with the
% three lowest points.
zp_tent=zp_sorted(1:3);
xp_tent=xp_sorted(1:3);
yp_tent=yp_sorted(1:3);

while ~isempty(zp_out)
    
    % find all points not in the tent that are within the tent's convex hull
    
    hullindices=convhull(xp_tent,yp_tent);
    out_within_tent=inpolygon(xp_out,yp_out,...
        xp_tent(hullindices),yp_tent(hullindices));
    
    % eliminate those points from future consideration
    
    zp_out(out_within_tent)=[];
    xp_out(out_within_tent)=[];
    yp_out(out_within_tent)=[];
    
    
    if ~isempty(zp_out);
        
        % Add the lowest point in the 'out' set. to the tent.
        zp_tent=[zp_tent zp_out(1)];
        xp_tent=[xp_tent xp_out(1)];
        yp_tent=[yp_tent yp_out(1)];
        
        % remove it from the 'out' set. 
        zp_out(1)=[];
        xp_out(1)=[];
        yp_out(1)=[];
    end
    
    % keep going until all points are in tent or excluded
end

%Store the bottom tent and the Delaunay triangulation of it for future
%display. 

%Add "spots" on the boundary
xp_tent = [xp_tent x1'];
yp_tent = [yp_tent y1'];
zp_tent = [zp_tent r'];

bottomtent=[xp_tent' yp_tent' zp_tent'];
bottomdt=DelaunayTri(xp_tent',yp_tent');



% the volume is calculated by building triscatteredinterp functions for
% both the top and the bottom and creating pixel by pixel maps of the top 
% and bottom tents.

% Using the 'natural' method because it gives a much smoother
% interpolation.
if numel(topdt.X(:,1))~=numel(toptent(:,1))
    Ftop=TriScatteredInterp(topdt,toptent(1:end-1,3),'natural');
    Fbottom=TriScatteredInterp(bottomdt,bottomtent(1:end-1,3),'natural');
else
    Ftop=TriScatteredInterp(topdt,toptent(:,3),'natural');
    Fbottom=TriScatteredInterp(bottomdt,bottomtent(:,3),'natural');
end

Xs=repmat(1:size(objmask,2),[size(objmask,1) 1]);
Ys=repmat((1:size(objmask,1))',[1 size(objmask,2)]);


% interpolate the tent heights at every pixel of the cell.
celltop=Ftop(Xs,Ys);
cellbottom=Fbottom(Xs,Ys);
% the interpolation procedure produces NaN for points outside the
% triangulation, so we set those to zero instead:
celltop(isnan(celltop))=0;
cellbottom(isnan(cellbottom))=0;

cellheight=celltop-cellbottom;
cellheight(~(objmask))=0;

cellvolume=sum(cellheight(:));