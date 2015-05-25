function [outpdf,x] = binomial_diff_pdf(p1,n1,p2,n2)

% This returns the PDF corresponding the random variable C = A-B, where A
% and B are random variables with a binomial PDF, with probabilities of p1,
% p2 and sample sizes of n1, n2 respectively.  The output variable x stores
% the values of the random variable C for convenience.

% First, let's make a PDF corresponding to A, with p1, n1.
pdf1 = binopdf(-n2:n1,n1,p1);

% Now let's get the PDF for B, with p2, n2.
pdf2 = binopdf(-n1:n2,n2,p2);


% Now let's get the PDF for -B.
pdf2 = fliplr(pdf2);  % First flip the pdf.  This goes from -n2 to n1

% Convolve to get the PDF of the sum
outpdf = conv(pdf1,pdf2);

% Domain is -2n1:2n2

x = (-2*n2):(2*n1);
