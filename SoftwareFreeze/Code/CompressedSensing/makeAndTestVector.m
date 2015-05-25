function [compOutVectors] = makeAndTestVector(M,vLength,nToChange)

v = rand(vLength,1)>.5;

for i = 1:nToChange
    v(floor(vLength/nToChange)*i) = 0;
end

Mv0 = M*v;

for i = 1:nToChange
    v(floor(vLength/nToChange)*i) = 1;
end

Mv1 = M*v;

compOutVectors = horzcat(Mv0,Mv1);