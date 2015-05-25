abund100 = load('abund100_from.01.txt');
mygenes = load('mygenes_ALLDATA_MatlabFriendly.txt');
%mygenes = sortrows(mygenes,2);

txt = ['gapdh' 'eef2'];

kList = [.01 .1 1 10 100 1000];
tauList = [.01 .1 .2 .3 .4 .5 .6 .7 .8 .9];

figure;
[hC hC] = contourf(abund100,20);
set(hC,'LineStyle','none');
xlabel('tau','fontSize',20);
ylabel('k','fontSize',20);
set(gca,'YTick',1:numel(kList))
set(gca,'XTick',1:numel(tauList))
set(gca,'YTickLabel',kList)
set(gca,'XTickLabel',tauList)

tmp = mygenes(:,2);
tmp = tmp*(log10(24)/6);
tmp(tmp<1) = tmp(tmp<1)+.5;
kData = tmp;

tmp = mygenes(:,5);
tmp = tmp*10;
tmp(tmp<1) = tmp(tmp<1)+1;
tauData = tmp;

hold on;
%scatter(mygenes(:,5)*10,mygenes(:,2)*(log10(24)/6),100,mygenes(:,6),'filled')
%scatter(mygenes(:,5)*10,mygenes(:,2)*(log10(24)/6),100,'k')

scatter(tauData,kData,100,mygenes(:,6),'filled')
scatter(tauData,kData,100,'k')

text(tauData(1)+.1,kData(1)+.1,'gapdh');
text(tauData(2)+.1,kData(2)+.1,'eef2');
text(tauData(3)+.1,kData(3)+.1,'actn4');
text(tauData(4)+.1,kData(4)+.1,'lmna');
text(tauData(5)+.1,kData(5)+.1,'tbcb');
text(tauData(6)+.1,kData(6)+.1,'supt5h');
text(tauData(7)+.1,kData(7)+.1,'icam1');
text(tauData(8)+.1,kData(8)+.1,'znf444');
text(tauData(9)+.1,kData(9)+.1,'slc1a5');
text(tauData(10)+.1,kData(10)+.1,'usf2');
text(tauData(11)+.1,kData(11)+.1,'pabpc1');
text(tauData(12)+.1,kData(12)+.1,'lum');
text(tauData(13)+.1,kData(13)+.1,'ubc');
text(tauData(14)+.1,kData(14)+.1,'gas6');
text(tauData(15)+.1,kData(15)+.1,'gaa');
text(tauData(16)+.1,kData(16)+.1,'ebf1');
text(tauData(17)+.1,kData(17)-.1,'rbm3');

%%

x = load('genesnear100.txt');

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

hold on;
scatter(x(:,4)*10,x(:,1)*(log10(24)/6),100,x(:,5),'filled')
scatter(x(:,4)*10,x(:,1)*(log10(24)/6),100,'k')

%%

x = load('121130_ActD_0h_UBC.txt');
%imagesc(flipud(expt'));
tmp = x(:,2)*(size(expt,1)/max(x(:,2)));

figure;
imagesc(expt);
hold on; scatter(x(:,3),tmp,100,'r','filled')
xlabel('mRNA','fontSize',20)
ylabel('Volume (scaled)','fontSize',20)
title('mu = mu0*vol','fontSize',20)