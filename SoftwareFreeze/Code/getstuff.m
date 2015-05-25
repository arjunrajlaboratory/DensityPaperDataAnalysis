function [] = getstuff(filename)

x = load(filename); 
filename
m=mean(x(:,4)); c=corr2(x(:,3),x(:,4)); [neg]=Intron_Spot_Frequency(x);
m
c
neg