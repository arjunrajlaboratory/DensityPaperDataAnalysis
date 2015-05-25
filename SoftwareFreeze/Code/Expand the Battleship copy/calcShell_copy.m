function [shell shell_filt] = calcShell(obj,x,y,z,expFactor,filtersize)

mask = imresize(obj.object_mask.mask,expFactor);

x2 = x*expFactor;
y2 = y*expFactor;
z2 = z*expFactor;

F = TriScatteredInterp(x2,y2,z2);

ti = 1:size(mask,1);
tf = 1:size(mask,2);
[qx,qy] = meshgrid(tf,ti);
shell = F(qx,qy);
shell_filt = medfilt2(shell,[filtersize filtersize]);