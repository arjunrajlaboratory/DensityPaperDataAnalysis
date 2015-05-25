% abund = 1000;
% kList = [1 5 10 50 100 500 1000];
% tauList = [.1 .2 .3 .4 .5 .6 .7 .8 .9 1];
 
abund = 100;
kList = [.01 .1 1 10 100 1000];
tauList = [.01 .1 .2 .3 .4 .5 .6 .7 .8 .9];

%abund = 788;
%kList = [.1 1 10];
%tauList = [ .7 .8 .9];

%abund = 10;
%kList = [.01 .1 1 10 100 1000];
%tauList = [.1 .2 .3 .4 .5 .6 .7 .8 .9 1];

clear rValMat

for i = 1:numel(tauList) %tau
    text = sprintf('on i=%d of %d',i,numel(tauList));
    disp(text);
    tau0 = tauList(i);
    for j = 1:numel(kList) %k
        text = sprintf('on j=%d of %d',j,numel(kList));
        disp(text);
        k0 = kList(j);
        
        mu0 = abund/tau0;
        lam0 = tau0*k0;
        gam0 = k0*(1-tau0);
        
        vArray = 0.3:.03:1.7;
        expt = 0;
        for p = 1:numel(vArray)
            disp(p)
            vol = vArray(p);
            
            %mu = m0*vol;
            %del = del0/vol;
            %lam = vol*gam0/(mu0-vol);
            %gam = lam0*(mu0/vol-1);
            
            cell = getpdf(lam0,gam0,mu0*vol,1);
            
            if expt == 0
                expt = cell;
            else
                expt = appendCell(expt,cell);
            end
        end
        [corr avgN] = calcCorr(expt);
        rValMat(j,i) = corr;
    end
end

%filename = sprintf('muVar_abund%d.mat',abund);
%save(filename,'rValMat');
%save('~/Dropbox/delMat.mat','rValMat');

figure;
%[x1 y1] = meshgrid(1:numel(tauList),1:numel(kList));
%scatter(x1(:),y1(:),50,rValMat(:),'filled');
[hC hC] = contourf(rValMat,20);
set(hC,'LineStyle','none');
xlabel('tau','fontSize',20);
ylabel('k','fontSize',20);
set(gca,'YTick',1:numel(kList))
set(gca,'XTick',1:numel(tauList))
set(gca,'YTickLabel',kList)
set(gca,'XTickLabel',tauList)

title('mu = mu0*vol','fontSize',20);