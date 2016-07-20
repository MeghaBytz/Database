function etchRate = calcEtchRate(radicals, ions, rxnProb, sputteringYield, ionStimulatedDesorption, wallLossFactors,radicalVelocities, ionVelocities)
%Surface Kinetics model
%determine free surface sites


%default to values of MgO but need to make this generalized
density = 3.56e+3; %kg/m^3
M = 40.3;
Na = 6.02e+23;


numberOfIons = length(ions);
numberOfRadicals = length(radicals);
chemicalEtchTerm = 0;
desorption = 0;
chemicalEtch = 0;
sputtering = 0;
for i = 1:numberOfRadicals
    radicalFlux(i) = radicals(i)*wallLossFactors(1)*radicalVelocities(i);
end

for i = 1:numberOfIons
    ionFlux(i) =  ions(i)*wallLossFactors(2)*ionVelocities(i);
end
for i = 1: numberOfRadicals
    chemicalEtchTerm = rxnProb(i)*radicalFlux(i) + chemicalEtchTerm;
end

for i = 1: numberOfIons
    desorption = ionStimulatedDesorption(i)*ionFlux(i) + desorption;
end

theta = chemicalEtchTerm/(desorption+chemicalEtchTerm);

for i = 1:numberOfRadicals
    (1-theta)*rxnProb(i)*radicalFlux(i) + chemicalEtch;
end

for i = 1:numberOfRadicals
    sputtering = (1-theta)*sputteringYield(i)*ionFlux(i) + sputtering;
end

etchRate = chemicalEtch + sputtering;
etchRate = etchRate*6e+11*M/(Na*density); %nm/min
%add unit conversions and velocity terms to make -- fluxes in units of 

end
    


