%%

obj = objects(1);
objmask = obj.object_mask.mask;

x = obj.channels.nir.spotCoordinates(:,2);
y = obj.channels.nir.spotCoordinates(:,1);
z = obj.channels.nir.spotCoordinates(:,3);

figure;
scatter3(x,y,z,'filled','m','MarkerEdgeColor','k');
hold on;

for i=1:20:900

idx = i;
%num = 1;

topline = celltop(idx,1:size(celltop,2));
sidxtop = min(find(topline~=0));
eidxtop = max(find(topline~=0));
topline(find(topline==0)) = [];

bottomline = cellbottom(idx,1:size(cellbottom,2));
sidxbottom = min(find(bottomline~=0));
eidxbottom = max(find(bottomline~=0));
bottomline(find(bottomline==0)) = [];
%figure; 
%line(sidxtop:eidxtop, topline,'LineWidth',2,'Color','k','LineSmoothing','on'); hold on;

plot3(repmat(idx,1,numel(sidxtop:eidxtop)),sidxtop:eidxtop,topline,'LineWidth',2,'Color','k','LineSmoothing','on');
plot3(repmat(idx,1,numel(sidxbottom:eidxbottom)),sidxbottom:eidxbottom,bottomline,'LineWidth',2,'Color','k','LineSmoothing','on');


b = bwboundaries(obj.object_mask.mask);
m = b{1};
plot(m(:,1),m(:,2),'LineWidth',2,'Color','b','LineSmoothing','on')
end
