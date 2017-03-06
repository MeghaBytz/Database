function [ionFlux radicalFlux] = calcFlux(ionDensities,radicalDensities,Te)

%inputs = density, electron Te,Ti distribution? 
%output = flux

global noUnknowns
global numberOfIons
global numberOfRadicals
global massRadicals
global massIons
global MO
global MO2
global Ti
%assumptions
%Radicals follow Maxwellian distribution
%Vertical component of ions follows bohm distribution. Azimuthal component
%is Maxwellian

massRadicals = MO;
massIons = [MO2 MO];
me=9.1095E-31;
MO=1836*16*me; % mass of an Oxygen atom
MO2=2*MO; % mass of an Oxygen Molecule
ee=1.6022E-19;

for i=1:numberOfIons
    ionFluxVz(i) = ionDensities(i)*sqrt((ee*Te/(massIons(i)))); %bohm velocity flux calculation 
    ionFluxVx(i) = ionDensities(i)*sqrt(8*ee*Ti*pi/(2*massIons(i)))/2; %cm^2/s
end
%calculate ion fluxes
for i=1:numberOfRadicals
    radicalFluxVz(i) = radicalDensities(i)*sqrt((8*ee*Ti/(pi*massRadicals(i))))/4;
    radicalFluxVx(i) = radicalDensities(i)*sqrt(8*ee*Ti*pi/(2*massRadicals(i)))/2; %cm^2/s
end

ionFlux = [ionFluxVz;ionFluxVx]*1e-4; %cm^2/s
radicalFlux = [radicalFluxVz;radicalFluxVx]*1e-4; %cm^2/s

