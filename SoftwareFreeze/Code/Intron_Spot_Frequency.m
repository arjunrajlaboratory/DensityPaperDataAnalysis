function [neg] = Intron_Spot_Frequency(x)

%x = load('Numbers.txt');

introns = x(:,6);
cyclin = x(:,5);

%cyclin pos

introns = introns(cyclin>100);  
introns = introns(introns<5);   
numint = numel(introns);

pos = mean(introns./4);

%cyclin neg

introns = x(:,6);
cyclin = x(:,5);

introns = introns(cyclin<100);  
introns = introns(introns<3);   
numint = numel(introns);

neg = mean(introns./2);