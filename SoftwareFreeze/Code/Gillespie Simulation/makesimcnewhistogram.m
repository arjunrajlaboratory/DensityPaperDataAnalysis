
% 8/14/05
% This file is to create a Gillespie simulation based on given reactions.
% It will create an .m file which will run the simulation.
% Additionally, it will create a function .m file
% which can be used with one of Matlab's ODE solvers
% to simulate the deterministic case.

% In general, the code will run fastest if the most frequent
% reactions are placed first, allowing for a quick selection of
% which reaction occurred.

function uniquetok = makesimcnewhistogram(inputfile, outputfile)

% First thing is to read in inputfile
fid = fopen(inputfile);

finc = fopen('basicgillespiehistogram.c');
foutc = fopen([outputfile 'histomex.c'],'w');

% The first line is the number of Gillespie steps to do
gillstepsline = fgetl(fid);
maxgillespiesteps = str2num(gillstepsline);

clear inlines lines rate ratenames ratelines

i=1;
while 1
  tline = fgetl(fid);  %Get line
  while isempty(tline) %ignore blank lines
    tline = fgetl(fid);
  end;
  
  if ~ischar(tline), break, end
  if strcmp(tline,'Initial values'), break, end
  inlines{i}=tline;
  i=i+1;
end;

i=1;
while 1  %This loop is to read all the initial values
  tline = fgetl(fid);
  while isempty(tline) %ignore blank lines
    tline = fgetl(fid);
  end;
  
  if ~ischar(tline), break, end
  initiallines{i} = tline;
  i=i+1;
end;

% Okay, now let's extract the rates from the beginning of the lines
j=1;
for i = 1:length(inlines)
  clear R K
  [ratelin,lin] = strread(inlines{i},'%s%s','delimiter',':');
  ratelines{i} = strrep(ratelin{1},',',' ');  %replace comma with space
  lines{j} = lin{1};
  [R,K] = strread(ratelines{i},'%s%f','delimiter','=');
  ratenames{j} = R{1};
  rate(j) = K(1);
  
  if length(R)==2  %If it's bidirectional...
    % First, split the line into the lhs and the rhs:
    j = j+1;  %Now we add the rates and the reverse rxn:
    ratenames{j} = R{2};
    rate(j) = K(2);
    [LHS,RHS] = strread(lin{1},'%s%s','delimiter','=');
    lines{j} = [RHS{1},'=',LHS{1}];
  end;
  j = j+1;
  
end;

fclose(fid);

fid = fopen([outputfile,'histo.m'],'w');

fprintf(fid,'%% First, let''s set all the simulation parameters\n');
fprintf(fid,'%sparams;\n\n',outputfile);

numrxns = length(lines);

currrxn = 1;

clear alllhs allrhs alltok rxntok

fprintf(fid,'for i = 1:maxgillespiesteps\n');

for currrxn = 1:numrxns
  
  clear lhstok rhstok;
  
  % First, split the line into the lhs and the rhs:
  [LHS,RHS] = strread(lines{currrxn},'%s%s','delimiter','=');
  
  % Now get all tokens from the LHS...
  if ~isempty(LHS)
    lhstok = strtrim(strread(LHS{1},'%s','delimiter','+'));
    rxntok{currrxn}.lhs = lhstok; %remove whitespace
  else
    lhstok = {};
    rxntok{currrxn}.lhs = {};
  end;
  
  
  % Now get all tokens from the LHS...
  if ~isempty(RHS)
    rhstok = strtrim(strread(RHS{1},'%s','delimiter','+'));
    rxntok{currrxn}.rhs = rhstok; 
    %remove whitespace
  else
    rhstok = {};
    rxntok{currrxn}.rhs = {};
  end;

  alltok{currrxn} = [lhstok',rhstok'];
end;

% Make a list of all the tokens we have in the whole thing
tok = [alltok{:}];

% This will just keep all unique tokens, yielding
% a list of species:
uniquetok = unique(tok);

% This gives the total number of species:
nspecies = length(uniquetok);

% Okay, now let's write the number of species and 
% number of reactions into the C code

lin = 'junk';
while ~strcmp(lin,'//INCLUDE NSPECIES HERE')
  lin = fgetl(finc);
  fprintf(foutc,'%s\n',lin);
end;
fprintf(foutc,'#define NSPECIES %d\n',nspecies);

while ~strcmp(lin,'//INCLUDE NUMRXNS HERE')
  lin = fgetl(finc);
  fprintf(foutc,'%s\n',lin);
end;
fprintf(foutc,'#define NUMRXNS %d\n',numrxns);


while ~strcmp(lin,'//INSERT ALL VARIABLE DECLARATIONS HERE')
  lin = fgetl(finc);
  fprintf(foutc,'%s\n',lin);
end;

fprintf(foutc,'  long %s',uniquetok{1});
for i = 2:nspecies
  fprintf(foutc,',%s',uniquetok{i});
end;
fprintf(foutc,';\n');

fprintf(foutc,'  double %s',ratenames{1});
for i = 2:length(ratenames)
  fprintf(foutc,',%s',ratenames{i});
end;
fprintf(foutc,';\n');



while ~strcmp(lin,'//UNPACK ALL SPECIES HERE')
  lin = fgetl(finc);
  fprintf(foutc,'%s\n',lin);
end;

for i = 1:nspecies
  fprintf(foutc,'  %s = (long)species[%d];\n',uniquetok{i},i-1);
end;



while ~strcmp(lin,'//UNPACK ALL RATES HERE')
  lin = fgetl(finc);
  fprintf(foutc,'%s\n',lin);
end;

for i = 1:length(ratenames)
  fprintf(foutc,'  %s = rates[%d];\n',ratenames{i},i-1);
end;


% Let's generate a list which, for each species, gives the
% LHSs it is involved in.  This is used later to determine
% which propensities to update.
clear lhsrxnlist

for i = 1:nspecies
  lhsrxnlist{i} = [];
  for currrxn = 1:numrxns
    if sum( ismember(rxntok{currrxn}.lhs, uniquetok(i)) );
      lhsrxnlist{i} = [lhsrxnlist{i},currrxn];
    end;
  end;
end;

% While we're at it, let's make a stoichiometry matrix.
% (Could be integrated into above loops, but I'll repeat for
% clarity...)
lhsstoichiometry = zeros(nspecies,numrxns);
rhsstoichiometry = zeros(nspecies,numrxns);
for i = 1:nspecies
  for j = 1:numrxns
    lhsstoichiometry(i,j) = sum( ismember(rxntok{j}.lhs, uniquetok(i)) );
    rhsstoichiometry(i,j) = sum( ismember(rxntok{j}.rhs, uniquetok(i)) );
  end;
end;

    
%***************** Okay, now we can really make some code************

% First, let's generate the code for timestep sampling:

%fprintf(fid,'alpha = 0\n');
fprintf(fid,'cumulativeprop = cumsum(propensity);\n');
fprintf(fid,'alpha = cumulativeprop(numrxns);\n');

fprintf(fid,'deltaT = -1/alpha*log(rand);\n');
fprintf(fid,'currT = currT + deltaT;\n');


% Now let's find out what actually happened:

fprintf(fid,'p = rand*alpha;\n');


% Save data

while ~strcmp(lin,'//INSERT SAVE HERE')
  lin = fgetl(finc);
  fprintf(foutc,'%s\n',lin);
end;



% Now let's generate code for saving some data:
for j=1:nspecies
  fprintf(fid,'savespecies(i,%d) = %s;\n',j,uniquetok{j});
  fprintf(foutc,'species_out[savecount*NSPECIES+%d] = %s;\n',j-1,uniquetok{j});
end;




while ~strcmp(lin,'//INSERT IF STATEMENT HERE')
  lin = fgetl(finc);
  fprintf(foutc,'%s\n',lin);
end;


for i = 1:numrxns
  % first, write code to check the appropriate rxn happened
  if i == 1
    fprintf(fid,'if p<cumulativeprop(%d)\n',i);
    fprintf(foutc,'if (p<cumpropensities[%d]) {\n',i-1);
  else
    fprintf(fid,'elseif p<cumulativeprop(%d)\n',i);
    fprintf(foutc,'} else if (p<cumpropensities[%d]) {\n',i-1);
  end;
  
  fprintf(foutc,'  // rxn: %s\n',lines{i});
  fprintf(fid,  '  %% rxn: %s\n',lines{i});

  % okay, we're in the if block for this reaction
  % now we must update the reactants...
  allstoichiometry = rhsstoichiometry-lhsstoichiometry;
  idx = find(allstoichiometry(:,i)~=0);
  for j = idx'
    fprintf(foutc,'  %s=%s + %d;\n',uniquetok{j},uniquetok{j},allstoichiometry(j,i));
    fprintf(fid,'  %s=%s + %d;\n',uniquetok{j},uniquetok{j},allstoichiometry(j,i));
  end;
  
  % alright, now let's find out which reactions will require
  % updated propensities...
  
  % this line finds all the species which change:
  idx = find( (lhsstoichiometry(:,i)-rhsstoichiometry(:,i)) ~= 0 );
  updatedrxns = unique([lhsrxnlist{idx}]);  % combine the list
  
  fprintf(foutc,'\n');
  fprintf(fid,'\n');
  

  % okay, now let's update the appropriate propensities
  for currrxn = updatedrxns
    fprintf(foutc,'  //update propensity for %s\n',lines{currrxn});
    fprintf(foutc,'  propensities[%d] = %s',currrxn-1,ratenames{currrxn});
    
    fprintf(fid,'  %%update propensity for %s\n',lines{currrxn});
    fprintf(fid,'  propensity(%d) = %s',currrxn,ratenames{currrxn});
      
    idx = find(lhsstoichiometry(:,currrxn)>0);
    for j = idx'
      fprintf(foutc,'*%s',uniquetok{j});
      fprintf(fid,'*%s',uniquetok{j});
      for k = 1:lhsstoichiometry(j,currrxn)-1
	fprintf(foutc,'*(%s-%d)',uniquetok{j},k);
	fprintf(fid,  '*(%s-%d)',uniquetok{j},k);
      end;
    end;
    fprintf(foutc,';\n');
    fprintf(fid,';\n');
    
  end;
  

end;

fprintf(foutc,'}\n');
fprintf(fid,'end;\n');


while 1  % write out the rest of the file
  tline = fgetl(finc);
  if ~ischar(tline),break,end
  fprintf(foutc,'%s\n',tline);
end;


fprintf(fid,'times(i) = currT;\n');

fprintf(fid,'end;\n');  % end of main loop

if fid~=1
  fclose(fid);
end;

fclose(finc);
fclose(foutc);

fprintf('Compiling %smex.c... ',outputfile);
tic;
eval(['mex ' outputfile 'histomex.c']);
%unix(['mex ' outputfile 'histomex.c']);
t = toc;
fprintf('elapsed time = %g seconds\n',t);

%**********************************************%

% Okay, now let's generate the parameters file.
% This allows one to set all the rates and so forth.
% It is called at the beginning of the simulation, and
% also computes the propensities ahead of time.

%fid = 1;

fid = fopen([outputfile,'params.m'],'w');

fprintf(fid,'%% Parameter file\n\n\n');

fprintf(fid,'%% Simulation parameter values\n\n');

fprintf(fid,'%% Maximum number of Gillespie steps\n');
fprintf(fid,'maxgillespiesteps = %d;\n\n',maxgillespiesteps);
fprintf(fid,'%% Initial time\n');
fprintf(fid,'currT = 0;\n\n\n');

fprintf(fid,'%% Initial values for species\n\n');
for i = 1:nspecies
  %fprintf(fid,'%s = ;\n',uniquetok{i});
  fprintf(fid,'%s;\n',initiallines{i});
end;

fprintf(fid,'\n\n%% EVERYTHING BELOW IS BOOKKEEPING; DO NOT ALTER!\n');

fprintf(fid,'\n\n%% Reaction rates\n\n');
for i = 1:numrxns
  fprintf(fid,'%% Rxn: %s\n',lines{i});
  fprintf(fid,'%s = %g;\n',ratenames{i},rate(i));
%  fprintf(fid,'rate(%d) = ;\n\n',i);
end;

% Let's put in some of the basic parameters:
fprintf(fid,'numrxns = %d;\n',numrxns);
fprintf(fid,'nspecies = %d;\n',nspecies);


% okay, first let's put the species numbers into the vector
fprintf(fid,'species = zeros(nspecies,1);\n');
for i = 1:nspecies
  fprintf(fid,'species(%d) = %s;\n',i,uniquetok{i});
end;

% Now let's put the reaction rates into the vector
fprintf(fid,'rates = zeros(numrxns,1);\n');
for i = 1:numrxns
%  fprintf(fid,'rates(%d) = %g;\n',i,rate(i));
  fprintf(fid,'rates(%d) = %s;\n',i,ratenames{i});
end;


% Now let's copy to y0 for the ode solver
%fprintf(fid,'\n\ny0=species;\n\n');
fprintf(fid,'y0 = zeros(nspecies,1);\n');
for i = 1:nspecies
  fprintf(fid,'y0(%d) = %s;\n',i,uniquetok{i});
end;


% okay, now let's intialize the propensities

fprintf(fid,'%% Intialize the propensities...\n');
for currrxn = 1:numrxns
  fprintf(fid,'propensity(%d) = %s',currrxn,ratenames{currrxn});
  
  idx = find(lhsstoichiometry(:,currrxn)>0);
  for j = idx'
    fprintf(fid,'*%s',uniquetok{j});
    for k = 1:lhsstoichiometry(j,currrxn)-1
      fprintf(fid,'*(%s-%d)',uniquetok{j},k);
    end;
  end;
  fprintf(fid,';\n');
end;

fprintf(fid,'\n');

% Let's now preallocate the memory to store all the species
%fprintf(fid,'%% Preallocate memory\n\n');
%fprintf(fid,'times   = zeros(maxgillespiesteps,1);\n');
%fprintf(fid,'savespecies = zeros(maxgillespiesteps,%d);\n',nspecies);


fclose(fid);

%****************************************************

% Okay, let's make an ode program as well...

s = [outputfile,'ode'];

fid = fopen([s,'.m'],'w');

fprintf(fid,'function dy = %s(t,y)\n\n',s);

fprintf(fid,'dy = zeros(%d,1);\n',nspecies);

for i = 1:nspecies
  %first, let's find all rxns for which 
  %this species is in the LHS/RHS:
  lhsidx = find(lhsstoichiometry(i,:)>0);
  rhsidx = find(rhsstoichiometry(i,:)>0);

  fprintf(fid,'dy(%d) = ',i);
  
  for currrxn = lhsidx
    %fprintf(fid,'-rate(%d)',currrxn);    
    fprintf(fid,'-%d*%g',lhsstoichiometry(i,currrxn),rate(currrxn));
    idx = find(lhsstoichiometry(:,currrxn)>0);
    for j = idx'
      for k = 1:lhsstoichiometry(j,currrxn)
	fprintf(fid,'*y(%d)',j);
      end;
    end;
  end;

  for currrxn = rhsidx
    %fprintf(fid,' + rate(%d)',currrxn);    
    fprintf(fid,' + %d*%g',rhsstoichiometry(i,currrxn),rate(currrxn));
    idx = find(lhsstoichiometry(:,currrxn)>0);
    for j = idx'
      for k = 1:lhsstoichiometry(j,currrxn)
	fprintf(fid,'*y(%d)',j);
      end;
    end;
  end;
  
  fprintf(fid,';\n');
  
end;

