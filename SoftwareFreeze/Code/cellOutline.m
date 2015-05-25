%%

load data010;
obj = objects(3);
color = 'nir';
expFactor = 1;
filtersize = 30;

mask = imresize(obj.object_mask.mask,expFactor);
dapiMask = obj.channels.dapi.processor.mask;
b = bwboundaries(dapiMask);

x = obj.channels.(color).spotCoordinates(:,2);
y = obj.channels.(color).spotCoordinates(:,1);
z = obj.channels.(color).spotCoordinates(:,3);

[volume x2t y2t z2t x2b y2b z2b celltop cellbottom goodpts_top goodpts_bottom] = findVol_Expand_copy(obj, color, x, y, z, expFactor);

ti = 1:5:size(mask,1);
tf = 1:5:size(mask,2);
[qx,qy] = meshgrid(tf,ti);
cthalf = celltop(1:5:end,1:5:end);
cbhalf = cellbottom(1:10:end,1:10:end);

xt = x(goodpts_top);
yt = y(goodpts_top);
zt = z(goodpts_top);

figure;
mesh(qx,qy,cthalf); hold on;
%mesh(qx,qy,cbhalf);
scatter3(xt,yt,zt,50,zt,'filled');
axis off
%line(b{1}(:,2),b{1}(:,1),'LineWidth',5,'Color','r')

%%

height = celltop - cellbottom;
height(isnan(height)) = 0;
height(~mask) = 0;
volumeWithNuc = sum(height(:));
height(dapiMask) = 0;
volume = sum(height(:))

top_filt = medfilt2(celltop,[filtersize filtersize]);
bottom_filt = medfilt2(cellbottom,[filtersize filtersize]);
cthalf = top_filt(1:10:end,1:10:end);

figure;
mesh(qx,qy,cthalf); hold on;
scatter3(xt,yt,zt,50,zt,'filled');

height = top_filt - bottom_filt;
height(isnan(height)) = 0;
height(~mask) = 0;
volumeWithNuc = sum(height(:));
height(dapiMask) = 0;
volume_filt = sum(height(:))

%%%%%%

xn = x*expFactor;
yn = y*expFactor;
zn = z*expFactor;

ti = 1:10:size(mask,1);
tf = 1:10:size(mask,2);
[qx,qy] = meshgrid(tf,ti);
cthalf = ctf(1:10:end,1:10:end);

figure;
mesh(qx,qy,cthalf); hold on;
scatter3(xn,yn,zn,50,zn,'filled');
