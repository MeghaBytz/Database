% starting pressures in mTorr (180W)
ppvec=[6.02 13.06 21.81 30.60 34.86 43.64 52.09 60.75]
Qvec=[4 8 20] % O2 flow rate in sccm
ccpPower = [ 100 200 500 800]
gammaOvec=[0.5 0.43 0.33 0.27 0.23 0.2 0.15 0.13]
ind = 1;
for i = 1:length(ppvec)
    for j = 1:length(Qvec)
        for r = 1:length(ccpPower)
            allSynExpConditions(ind,1) = ppvec(i);
            allSynExpConditions(ind,2) = Qvec(j);
            allSynExpConditions(ind,3) = gammaOvec(i);
            allSynExpConditions(ind,4) = ccpPower(r);
            ind = ind+1;
        end
    end
end

%make synthetic etch data
for k = 1:length(allSynExpConditions)
    allSynData(k) = allSynExpConditions(k,4)*10/sqrt(allSynExpConditions(k,1));
    if allSynExpConditions(i,2)>12
        allSynData(k) = allSynData(k)/allSynExpConditions(i,2);
    end
end
