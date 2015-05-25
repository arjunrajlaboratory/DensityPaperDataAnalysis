y=load('densitymatrix_withNuclearVolume.txt');
m=find(y(:,5)==0);
y(m,:)=[];
figure; scatter(y(:,3),y(:,5));

x = load('densitymatrix.txt');
m=find(x(:,5)==0);
x(m,:)=[];
scatter(x(:,3),x(:,5))

y=load('densitymatrix_noBoundaryEnforce_includeNucVolume.txt');
m=find(y(:,5)==0);
y(m,:)=[];
figure; scatter(y(:,3),y(:,5));

%polyfit

%%%%%%%%%%%%%

x=load('densitymatrix_withNuclearVolume.txt');
m=find(x(:,5)==0);
x(m,:)=[];
corr2(y(:,3),y(:,5))

y=load('densitymatrix_withNuclearVolume.txt');
m=find(y(:,5)==0);
y(m,:)=[];
corr2(y(:,3),y(:,5))

m=find(y(:,3)>6E5);
y(m,:)=[];
p1=polyfit(y(:,3),y(:,5),1)
corr2(y(:,3),y(:,5))

y=load('densitymatrix_withNuclearVolume.txt');
m=find(y(:,5)==0);
y(m,:)=[];
m=find(y(:,3)<6E5);
y(m,:)=[];
p2=polyfit(y(:,3),y(:,5),1)
corr2(y(:,3),y(:,5))