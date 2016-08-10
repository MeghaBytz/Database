% starting pressures in mTorr (180W)
ppvec=[6.02 13.06 21.81 30.60 34.86 43.64 52.09 60.75]
Qvec=[40 50 50 50 50 50 50 50] % O2 flow rate in sccm
gammaOvec=[0.5 0.43 0.33 0.27 0.23 0.2 0.15 0.13]
ind = 1;
for i = 1:length(ppvec)
    for j = 1:length(Qvec)
        expConditions(ind,1) = ppvec(i);
        expConditions(ind,2) = Qvec(j);
        expConditions(ind,3) = gammaOvec(i);
        ind = ind+1;
    end
end