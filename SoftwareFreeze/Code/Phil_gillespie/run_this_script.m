% Run this script to generate the mex code.
% All times are in seconds.



% The reaction scheme is kept in mrnaburst.txt
% This command generates a bunch of files that can be used to run a Gillespie simulation.
tok = makesimcnewhistogram('mrnaburst.txt','burst');
% The variable tok contains the various species.  The order of these species is important,
% since the output vector has the species kept in this order.

% One of the files generated is burstparams.m.  It initializes various parameters for 
% the simulation itself.
burstparams

% The main code itself is written into a mex file called bursthistomex.
% This runs the Gillespie simulation itself, using the variables:
% species, rates, propensity
% all of which are defined in burstparams.  These are the initial conditions, the rates,
% and the initial propensities.  The program itself takes another couple parameters.
% The first one is sort of a dummy parameter of the initial time (not very useful).
% Another is the random number seed (I use sum(clock*100)).  The next two parameters are
% the number of data points to save (equally spaced in time) and the stop time.
% For example, if I want to save my data every 5 seconds for 5000 seconds total run time,
% I would do:


[t,s] = bursthistomex(0,species,rates,propensity,sum(clock*100),1000,10000);

% To plot, just type

plot(t,s);
legend(tok);
mrnaAtEnd = s(4,end)

% The reason I don't save every single iteration of the Gillespie algorithm is that this
% can rapidly saturate memory, and often you only need "snapshots" at regular intervals
% (for making histograms, for instance).  You can also try to run the generated file
% bursthisto.m, which contains a MATLAB script that runs the Gillespie algorithm but 
% saves data at every iteration (runs up to maxgillespiesteps number of iterations).  I 
% haven't used that in a long time, so I can't guarantee that it works, but you can try it.
% Finally, the program also spits out burstode.m, which is an ODE function you can solve
% using any of MATLAB's standard solvers like ode45.