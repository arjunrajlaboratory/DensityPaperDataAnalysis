%load data010;
%obj = objects(3);
obj = objects(1);
color = 'nir';
expFactor = 1;
filtersize = 30;

mask = imresize(obj.object_mask.mask,expFactor);
dapiMask = obj.channels.dapi.processor.mask;
b = bwboundaries(dapiMask);

x = obj.channels.(color).spotCoordinates(:,2);
y = obj.channels.(color).spotCoordinates(:,1);
z = obj.channels.(color).spotCoordinates(:,3);

[xt yt zt xb yb zb celltop cellbottom] = findVol_Expand2_getOutline(obj, color, x, y, z, expFactor,35);

celltop2 = medfilt2(celltop,[20,20]);

ti = 1:5:size(mask,1);
tf = 1:5:size(mask,2);
[qx,qy] = meshgrid(tf,ti);
cthalf = celltop2(1:5:end,1:5:end);
cbhalf = cellbottom(1:10:end,1:10:end);

%xt = x(goodpts_top);
%yt = y(goodpts_top);
%zt = z(goodpts_top);

figure;
mesh(qx,qy,cthalf); %hold on;
%surf(qx,qy,cthalf); hold on;
%mesh(qx,qy,cbhalf);
%scatter3(xt,yt,zt,50,zt,'filled');
axis off
%line(b{1}(:,2),b{1}(:,1),'LineWidth',5,'Color','r')