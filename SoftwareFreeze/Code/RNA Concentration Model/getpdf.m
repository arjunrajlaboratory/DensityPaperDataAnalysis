function pdf = getpdf(lam,gam,mu,del)


M = mastereqnmatrix(lam,gam,mu,del);

M(6,1) = 1;

sz = size(M);

b = zeros(sz(1),1);

b(6) = 1;
x = M\b;

x = x/sum(x);

x2 = reshape(x,2,length(x(:))/2);
pdf = sum(x2);

