
%function M = mastereqnmatrix(k1,ka,ki)
function M = mastereqnmatrix(lam,gam,mu,del)


k2 = 0;  %k2 is the transcription rate in the inactive state


B1 = [-gam lam; gam -lam];
B2 = diag([-mu -k2]);

%N = 1000*max([factor,1]);
N = round(4*max(mu,100));

M = spalloc(N,N,6*N);

i = 1;
M(i:i+1,i:i+1) = B1+B2-(i-1)/2;

for i = 3:2:N
  M(i:i+1,i:i+1) = B1+B2-del*eye(2)*(i-1)/2;
  M(i-2:i-1,i:i+1) = del*eye(2)*(i-1)/2;
  M(i:i+1,i-2:i-1) = -B2;
end;


