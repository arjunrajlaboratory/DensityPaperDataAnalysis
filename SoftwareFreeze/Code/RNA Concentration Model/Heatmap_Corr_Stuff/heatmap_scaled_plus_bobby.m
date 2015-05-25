abund100 = load('abund100_from.01.txt');
abund100 = abund100(3:end,:);
mygenes = load('mygenes_ALLDATA_MatlabFriendly.txt');
bobbygenes = load('bobbygenes_ALLDATA_MatlabFriendly.txt');

kList = [1 10 100 1000];
tauList = [.01 .1 .2 .3 .4 .5 .6 .7 .8 .9];

figure;
[hC hC] = contourf(abund100,20);
set(hC,'LineStyle','none');
xlabel('Duty cycle','fontSize',20);
ylabel('C/(degradation rate)','fontSize',20);
set(gca,'YTick',1:numel(kList))
set(gca,'XTick',1:numel(tauList))
set(gca,'YTickLabel',kList)
set(gca,'XTickLabel',tauList,'FontSize',15)

cb = colorbar;
set(cb,'FontSize',15);

tmp = mygenes(:,2);
tmp = log10(tmp)*(4/log10(24));
kData = tmp;

tmp = bobbygenes(:,1);
tmp = log10(tmp)*(4/log10(24));
kDataBobby = tmp;

tmp = mygenes(:,5);
tmp = tmp*10+1;
tauData = tmp;

tmp = bobbygenes(:,4);
tmp = tmp*10+1;
tauDataBobby = tmp;

hold on;
%scatter(mygenes(:,5)*10,mygenes(:,2)*(log10(24)/6),100,mygenes(:,6),'filled')
%scatter(mygenes(:,5)*10,mygenes(:,2)*(log10(24)/6),100,'k')

scatter(tauData,kData,200,mygenes(:,6),'filled')
scatter(tauData,kData,200,'k','LineWidth',1.2)
scatter(tauDataBobby,kDataBobby,200,'k','LineWidth',1.2)

text(tauData(1)+.1,kData(1)+.1,'GAPDH','FontSize',14,'LineWidth',2);
text(tauData(2)-.5,kData(2)-.1,'EEF2','FontSize',14);
text(tauData(3)+.1,kData(3)+.1,'ACTN4','FontSize',14);
text(tauData(4)+.1,kData(4)+.1,'LMNA','FontSize',14);
text(tauData(5)+.1,kData(5)+.1,'TBCB','FontSize',14);
text(tauData(6)+.1,kData(6)+.1,'SUPT5H','FontSize',14);
text(tauData(7)+.1,kData(7)+.1,'ICAM1','FontSize',14);
text(tauData(8)+.1,kData(8)+.1,'ZNF444','FontSize',14);
text(tauData(9)+.1,kData(9)+.1,'SLC1A5','FontSize',14);
text(tauData(10)+.1,kData(10)-.1,'USF2','FontSize',14);
text(tauData(11)-.7,kData(11)-.1,'PABPC1','FontSize',14);
text(tauData(12)+.1,kData(12)+.1,'LUM','FontSize',14);
text(tauData(13)+.1,kData(13)+.1,'UBC','FontSize',14);
text(tauData(14)+.1,kData(14)+.1,'GAS6','FontSize',14);
text(tauData(15)+.1,kData(15)+.1,'GAA','FontSize',14);
text(tauData(16)+.1,kData(16)+.1,'EBF1','FontSize',14);
text(tauData(17)+.1,kData(17)-.1,'RMB3','FontSize',14);

%colormap('cool'); brighten(-.6)