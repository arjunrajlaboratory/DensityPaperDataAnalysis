%lam = .01;
%mu = 100000;
%gam = 1;
%del = 1;

%abund = lam*mu/(lam+gam)
clear abundMat
param = [.01 .1 1 10 100 1000 100000 100000];

for i = 1:numel(param)
    lam = param(i);
    for j = 1:numel(param);
        mu = param(j);
        for k = 1:numel(param)
            gam = param(k);
            abundMat(i,j,k) = lam*mu/(lam+gam);
        end
    end
end

[x1 y1 z1] = meshgrid(1:numel(param),1:numel(param),1:numel(param));
scatter3(x1(:),y1(:),z1(:),50,abundMat(:),'filled');

m = 1;
for i = 1:numel(param)
    lam = param(i);
    for j = 1:numel(param)
        gam = param(j);
        tr(1) = lam/(lam+gam);
        m = m+1;
    end
end
tr = unique(tr');

for i = 1:numel(param)
    mu = param(i);
    for j = 1: numel(tr)
        trate = param(j);
        