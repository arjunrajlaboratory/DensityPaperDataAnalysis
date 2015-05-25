function [xf yf zf] = fillVol2(obj, color, ctf, cbf, nRealSpots, expFactor);

flag = 0;

% Resize cell top, bottom
ctf = imresize(ctf,expFactor);
ctf = ctf*expFactor;
cbf = imresize(cbf,expFactor);
cbf = cbf*expFactor;

[xf yf zf] = genUnifPtsExact(obj, color, ctf, cbf, nRealSpots, expFactor);

%if numel(xf) > nRealSpots*1.1
%    disp('too few spots');
%    return;
%end

%We want the number of spots to be within 5% of the real number of spots
%while abs(1-nRealSpots/numel(xf)) > 0.05
%    if numel(xf) > nRealSpots
%        flag = flag + 1;
%        if flag > 2
%            return;
%        end
%    end
%    [xf yf zf] = genUnifPts2(obj, color, ctf, cbf, nInitSpots, expFactor);
%    nInitSpots = nInitSpots + 1000;
%end