% Intron Frequency
x = load('NumbersVolume_GAPDHVol.txt');
cyclinPos = x(x(:,4)>100,:);
cyclinNeg = x(x(:,4)<100,:);

intPos = cyclinPos(:,5);
intPos(intPos>4)=[];
intNeg = cyclinNeg(:,5);
intNeg(intNeg>2)=[];

meanPos = mean(intPos)
errorPos = 1.96*std(intPos)/sqrt(numel(intPos))

meanNeg = mean(intNeg)
errorNeg = 1.96*std(intNeg)/sqrt(numel(intNeg))

% Intron Intensities
x = load('Vol_GAPDH_Intens_Cyclin.txt');
cyclinPos = x(x(:,5)>100,:);
cyclinNeg = x(x(:,5)<100,:);

intPos = cyclinPos(:,4);
intPos(intPos>4)=[];
intNeg = cyclinNeg(:,4);
intNeg(intNeg>2)=[];

meanPos = mean(intPos)
errorPos = 1.96*std(intPos)/sqrt(numel(intPos))

meanNeg = mean(intNeg)
errorNeg = 1.96*std(intNeg)/sqrt(numel(intNeg))