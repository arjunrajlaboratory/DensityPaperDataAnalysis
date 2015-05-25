function [ cellvolume cellheight celltop cellbottom] = tentconstruct_fromRNAspots_smallRegion( dapiMask, objmask , rnacoords , image )
% obj is the image object representing for example a cell. 
% rnacoords is an Mx3 with the [x y z] coordinates of every RNA spot to be
% used in constructing the top and bottom envelopes (tents) around the
% object. 



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
bottomtent=[xp_tent' yp_tent' zp_tent'];
bottomdt=DelaunayTri(xp_tent',yp_tent');



% the volume is calculated by building triscatteredinterp functions for
% both the top and the bottom and creating pixel by pixel maps of the top 
% and bottom tents.

% Using the 'natural' method because it gives a much smoother
% interpolation.
Ftop=TriScatteredInterp(topdt,toptent(:,3),'natural');
Fbottom=TriScatteredInterp(bottomdt,bottomtent(:,3),'natural');

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

%Exclude region around nucleus.
cellheight(dapiMask)=0;

cellvolume=sum(cellheight(:));

%end

% ----------------- Now plot! ----------------- %


% plotrows = 2;
% plotcols = 3;
% 
% 
% figure(1);clf;
% subplot(plotrows,plotcols,1)
% imshow(image)
% hold on
% 
% heightmap=jet(max(zp)-min(zp)+1);
% 
% for i=1:length(zp)
% 
% line(xp(i),yp(i),'Marker','.','MarkerEdgeColor',...
%     heightmap(zp(i)-min(zp)+1,:))
% end
% 
% hold off
% 
% figure(1);subplot(plotrows,plotcols,2)
% imshow(image)
% hold on
% dt=DelaunayTri(xp_tent',yp_tent');
% hold on
% triplot(dt)    
% 
% for i=1:length(zp_tent)
% line(xp_tent(i),yp_tent(i),'Marker','.','MarkerEdgeColor',...
%     heightmap(zp_tent(i)-min(zp)+1,:))
% end
% for i=1:length(zp_out)
% line(xp_out(i),yp_out(i),'Marker','.','MarkerEdgeColor',...
%     heightmap(zp_out(i)-min(zp)+1,:))
% end
% 
% hold off
% 
% figure(1);subplot(plotrows,plotcols,3)
% imshow(image)
% hold on
% dt=DelaunayTri(xp_tent',yp_tent');
% hold on
% triplot(dt)
% 
% for i=1:length(zp_tent)
% line(xp_tent(i),yp_tent(i),'Marker','.','MarkerEdgeColor',...
%     heightmap(zp_tent(i)-min(zp)+1,:))
% end
% for i=1:length(zp_out)
% line(xp_out(i),yp_out(i),'Marker','.','MarkerEdgeColor',...
%     heightmap(zp_out(i)-min(zp)+1,:))
% end
% 
% 
% % Need to now create topography maps. We can start with plotting the mask
% 
% figure(1);subplot(plotrows,plotcols,4)
% 
% imshow(objmask)
% 
% % Topography plot
% 
% figure(2);clf;subplot(1,plotcols,1)
% imshow((celltop-min(zp))/max(zp));colormap jet
% subplot(1,plotcols,2)
% imshow((celltop-min(zp))/max(zp));colormap jet
% 
% subplot(1,plotcols,3)
% heightZs=celltop-cellbottom;
% imshow((heightZs)/(max(zp)-min(zp)));colormap jet
% 
% 
% heightZs(~(objmask))=0;
% heightZs(isnan(heightZs))=0;
% heightZs=celltop-cellbottom;
% figure;
% imagesc((heightZs)/(max(zp)-min(zp)));colormap jet;