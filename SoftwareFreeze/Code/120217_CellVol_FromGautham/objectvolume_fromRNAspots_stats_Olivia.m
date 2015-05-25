function [ cellvolume_list ] = objectvolume_fromRNAspots_stats( obj, rna_color ,numRNA_subsample, numsamples)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

xp=obj.channels.(rna_color).fitdataRNAonly.yp_fit;
yp=obj.channels.(rna_color).fitdataRNAonly.xp_fit;
zp=obj.channels.(rna_color).fitdataRNAonly.rawzp;

numRNA=length(xp);

rnacoords=[xp', yp', zp'];

cellvolume_list=double(zeros([1 numsamples]));

for i=1:numsamples
    allindices=randperm(numRNA);
    
    if numRNA_subsample<numRNA
        someindices=allindices(1:numRNA_subsample);
    else
        someindices=allindices
    end
    
    rnacoords_subsample=rnacoords(someindices,:);
    
    cellvolume_list(i)=tentconstruct_fromRNAspots(obj.object_mask.mask,rnacoords_subsample);
    
end

