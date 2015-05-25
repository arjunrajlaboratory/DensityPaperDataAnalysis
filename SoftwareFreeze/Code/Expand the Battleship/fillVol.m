function [xf yf zf] = fillVol(obj, color, ctf, cbf, nRealSpots, expFactor);

flag = 0;
nInitSpots = 1000;
[xf yf zf] = genUnifPts(obj, color, ctf, cbf, nInitSpots, expFactor);

if numel(xf) > nRealSpots*1.1
    disp('too few spots');
    return;
end

%We want the number of spots to be within 10% of the real number of spots
while abs(nRealSpots-numel(xf)) > 0.1*nRealSpots
    if numel(xf) > nRealSpots
        flag = flag + 1;
        if flag > 2
            return;
        end
    end
    [xf yf zf] = genUnifPts(obj, color, ctf, cbf, nInitSpots, expFactor);
    nInitSpots = nInitSpots + 1000;
end