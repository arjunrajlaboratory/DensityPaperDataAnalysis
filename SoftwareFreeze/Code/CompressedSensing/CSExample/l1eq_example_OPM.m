% Compressed sensing example

% Test out l1eq code (l1 minimization with equality constraints).
%
% Written by: Justin Romberg, Caltech
% Email: jrom@acm.caltech.edu
% Created: October 2005
%

% To reproduce the example in the documentation, uncomment the 
% two lines below
%load RandomStates
%rand('state', rand_state);
%randn('state', randn_state);

%% First example, reproduce random spikes

% signal length
N = 512;
% number of spikes in the signal
T = 20;
% number of observations to make
K = 120;

% random +/- 1 signal
x = zeros(N,1);
q = randperm(N);
x(q(1:T)) = sign(randn(T,1));
figure; subplot(3,1,1); plot(x); title('original signal');

%%

% measurement matrix
disp('Creating measurment matrix...');
A = randn(K,N);
A = orth(A')';
disp('Done.');
	
% observations
y = A*x;

subplot(3,1,2); plot(y,'k'); title('measurements');

%%

% initial guess = min energy
x0 = A'*y;

% solve the LP
tic
xp = l1eq_pd(x0, A, [], y, 1e-3);
toc


subplot(3,1,3); plot(xp,'r'); title('reconstructed signal');

% large scale
% Afun = @(z) A*z;
% Atfun = @(z) A'*z;
% tic
% xp = l1eq_pd(x0, Afun, Atfun, y, 1e-3, 30, 1e-8, 200);
% toc

%%

norm(x-xp)

%% Example 2, sine waves

% Try a combination of sin waves.
% THIS SIGNAL IS NOT SPARSE!!! THAT'S WHY IT DOESN'T WORK.

x = 1/2*sin(pi*((1:512)/512)*10) + 1/2*sin(pi*((1:512)/512)*30);
x = x';
figure; subplot(3,1,1); plot(x); title('original signal');

%%
y = A*x;

subplot(3,1,2); plot(y,'k'); title('measurements');
sound(repmat(x,10,1))

%%


% initial guess = min energy
x0 = A'*y;

% solve the LP
tic
xp = l1eq_pd(x0, A, [], y, 1e-3);
toc


%plot(x)
%hold on;
%plot(xp,'r');

sound(repmat(x,20,1))
sound(repmat(xp,20,1))


subplot(3,1,3); plot(xp,'r'); title('reconstructed signal');

%%

%%%% Let's do it on the DCT of the function
% This turns the signal into a sparse signal.
x = 1/2*sin(pi*((1:512)/512)*10) + 1/2*sin(pi*((1:512)/512)*30);
x = x';

figure; subplot(4,1,1); plot(x); title('original signal');
sound(repmat(x,10,1))

%%

dtx = dct(x);
subplot(4,1,2); plot(dtx); title('transformed signal (direct cosine transform)');

%%

y = A*dtx;

x0 = A'*y;

xp = l1eq_pd(x0, A, [], y, 1e-3);

%plot(dtx);
%hold on;
%plot(xp,'r');
subplot(4,1,3); plot(xp,'r'); title('reconstructed transformed signal');

%%

x2 = idct(xp);

%hold off;
subplot(4,1,4);
plot(x)
hold on;
plot(x2,'r');
title('reconstructed original signal');
sound(repmat(x,20,1))
sound(repmat(x2,20,1))

