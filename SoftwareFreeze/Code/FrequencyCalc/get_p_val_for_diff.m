function p_val = get_p_val_for_diff(frequency1, n1, frequency2, n2)

% If you have frequency1 frequency of spots on n1 potential gene copies on the
% intact chromosome, and frequency2 frequency of spots on n2 potential gene
% copies on the fragment, this function returns the p-value for the
% difference in frequencies being meaningful.
% For example: I have 10 cells, and I see 0.3 frequency of spots on the 20
% potential intact copies.  On the other 10 gene copies on the fragment,
% I see a 0.8 frequency of spots.  To get the p-value for this difference,
% I type:
% get_p_val_for_diff(0.3,20,0.8,10)

% The null hypothesis is that both have the same frequency.  This frequency
% is given by the weighted average:

p = (frequency1*n1 + frequency2*n2)/(n1+n2);

p_val = binomial_diff_mean_test(p,n1,n2,abs(frequency1-frequency2));


