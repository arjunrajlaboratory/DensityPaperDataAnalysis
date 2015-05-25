function pdf = calc_burst_pdf(p1,p2);

%taken from supplement of plos bio paper

%parameter p1 = mu/gamma;
%parameter p2 = lambda/delta;

pdf = 0;

sz = round(4*max(p1,100));

for m = 1:sz
    
    x = ( gamma(p2+m)/(gamma(p2)*factorial(m)) ) * ( (1 - 1/(p1+1))^m ) * (1/(p1+1))^p2;
    if x == Inf | isnan(x)
        break;
    end
    
    pdf(m) = x;
    
end