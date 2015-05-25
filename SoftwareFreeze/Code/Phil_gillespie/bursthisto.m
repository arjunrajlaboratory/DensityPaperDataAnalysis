% First, let's set all the simulation parameters
burstparams;

for i = 1:maxgillespiesteps
cumulativeprop = cumsum(propensity);
alpha = cumulativeprop(numrxns);
deltaT = -1/alpha*log(rand);
currT = currT + deltaT;
p = rand*alpha;
savespecies(i,1) = P;
savespecies(i,2) = Pactive;
savespecies(i,3) = Pinactive;
savespecies(i,4) = mRNA;
if p<cumulativeprop(1)
  % rxn: Pactive=Pinactive
  Pactive=Pactive + -1;
  Pinactive=Pinactive + 1;

  %update propensity for Pactive=Pinactive
  propensity(1) = gam *Pactive;
  %update propensity for Pinactive=Pactive
  propensity(2) = lam*Pinactive;
  %update propensity for Pinactive=Pinactive+mRNA
  propensity(3) = k2 *Pinactive;
  %update propensity for Pactive=Pactive+mRNA
  propensity(4) = mu *Pactive;
elseif p<cumulativeprop(2)
  % rxn: Pinactive=Pactive
  Pactive=Pactive + 1;
  Pinactive=Pinactive + -1;

  %update propensity for Pactive=Pinactive
  propensity(1) = gam *Pactive;
  %update propensity for Pinactive=Pactive
  propensity(2) = lam*Pinactive;
  %update propensity for Pinactive=Pinactive+mRNA
  propensity(3) = k2 *Pinactive;
  %update propensity for Pactive=Pactive+mRNA
  propensity(4) = mu *Pactive;
elseif p<cumulativeprop(3)
  % rxn: Pinactive=Pinactive+mRNA
  mRNA=mRNA + 1;

  %update propensity for mRNA=
  propensity(5) = del *mRNA;
  %update propensity for mRNA=mRNA+P
  propensity(6) = mu_p *mRNA;
elseif p<cumulativeprop(4)
  % rxn: Pactive=Pactive+mRNA
  mRNA=mRNA + 1;

  %update propensity for mRNA=
  propensity(5) = del *mRNA;
  %update propensity for mRNA=mRNA+P
  propensity(6) = mu_p *mRNA;
elseif p<cumulativeprop(5)
  % rxn: mRNA=
  mRNA=mRNA + -1;

  %update propensity for mRNA=
  propensity(5) = del *mRNA;
  %update propensity for mRNA=mRNA+P
  propensity(6) = mu_p *mRNA;
elseif p<cumulativeprop(6)
  % rxn: mRNA=mRNA+P
  P=P + 1;

  %update propensity for P=
  propensity(7) = del_p *P;
elseif p<cumulativeprop(7)
  % rxn: P=
  P=P + -1;

  %update propensity for P=
  propensity(7) = del_p *P;
end;
times(i) = currT;
end;
