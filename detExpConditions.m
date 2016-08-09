% starting pressures in mTorr (180W)
ppvec=[6.02 13.06 21.81 30.60 34.86 43.64 52.09 60.75]
Qvec=[40 50 50 50 50 50 50 50] % O2 flow rate in sccm
ind = 1;
for i = 1:length(ppvec)
    for j = 1:length(Qvec)
        expConditions(ind,1) = ppvec(i);
        expConditions(ind,2) = Qvec(j);
        ind = ind+1;
    end
end