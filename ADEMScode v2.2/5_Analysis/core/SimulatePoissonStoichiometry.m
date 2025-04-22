%simulate Poisson stoichiometry data to test periodicity code.
function [stoichgt,stoichnoisy]=SimulatePoissonStoichiometry(periodgt,mu,noiseamp,listlen,plotflag)

stoichgt=[];
pd = makedist('Poisson',mu/periodgt);
for i=1:listlen
stoichgt=[stoichgt,random(pd)*periodgt];
end
stoichgt=stoichgt(stoichgt>0);  %truncated Poisson distribution

stoichnoisy=[];
for i=1:length(stoichgt)
stoichnoisy=[stoichnoisy,stoichgt(i)+random('Normal',0,noiseamp*sqrt(1+stoichgt(i)/12))];
end
stoichnoisy=stoichnoisy(stoichnoisy>0);

if plotflag
    figure;
    histogram(stoichgt);hold on
    histogram(stoichnoisy);
end
end