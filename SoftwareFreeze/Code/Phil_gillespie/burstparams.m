% Parameter file


% Simulation parameter values

% Maximum number of Gillespie steps
maxgillespiesteps = 5000;

% Initial time
currT = 0;


% Initial values for species

mRNA = 0;
Pactive = 1;
Pinactive = 0;
P = 100;


% EVERYTHING BELOW IS BOOKKEEPING; DO NOT ALTER!


% Reaction rates

% Rxn: Pactive=Pinactive
gam  = 0.005;
% Rxn: Pinactive=Pactive
lam = 0.002;
% Rxn: Pinactive=Pinactive+mRNA
k2  = 0;
% Rxn: Pactive=Pactive+mRNA
mu  = 0.03;
% Rxn: mRNA=
del  = 0.008;
% Rxn: mRNA=mRNA+P
mu_p  = 0.5;
% Rxn: P=
del_p  = 0.005;
numrxns = 7;
nspecies = 4;
species = zeros(nspecies,1);
species(1) = P;
species(2) = Pactive;
species(3) = Pinactive;
species(4) = mRNA;
rates = zeros(numrxns,1);
rates(1) = gam ;
rates(2) = lam;
rates(3) = k2 ;
rates(4) = mu ;
rates(5) = del ;
rates(6) = mu_p ;
rates(7) = del_p ;
y0 = zeros(nspecies,1);
y0(1) = P;
y0(2) = Pactive;
y0(3) = Pinactive;
y0(4) = mRNA;
% Intialize the propensities...
propensity(1) = gam *Pactive;
propensity(2) = lam*Pinactive;
propensity(3) = k2 *Pinactive;
propensity(4) = mu *Pactive;
propensity(5) = del *mRNA;
propensity(6) = mu_p *mRNA;
propensity(7) = del_p *P;

