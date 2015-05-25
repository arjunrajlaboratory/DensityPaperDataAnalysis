function p_val = binomial_diff_mean_test(p,n1,n2,dif)

% This function gives the p-value for rejecting the null-hypothesis of:
% if you toss n1+n2 coins with bias p, what are the chances of getting a
% difference in the mean between batch n1 with batch n2 that is greater
% than dif?

% Get the pdf of n1:
pdf1 = binopdf(0:n1,n1,p);
x1 = (0:n1)/n1;  % This is the range of values

% Same for n2:
pdf2 = binopdf(0:n2,n2,p);
x2 = (0:n2)/n2;

% Start with 
p_val = 0;
for i = 1:n1
    for j = 1:n2
        if abs(x1(i)-x2(j)) >= dif-0.001  % This extra 0.001 is to handle rounding errors
            p_val = p_val + pdf1(i)*pdf2(j);
        end;
    end;
end;
