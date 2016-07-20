% Integrate ODE's for global oxygen with new h_l factors
% Using the differential equations with average densities
% assumption of parabolic profiles for neutrals
% including O2, O2+, O+, O-, O & O2*
% Reaction coefficients have been corrected.

clear
global ee me MO MO2 ng Tg Ti l_p gammaO gammaO2m
global Efactor_O2 Efactor_O EnergyO2 sigO2 EnergyO sigO
global Krec Krec2 Krec3 Krec4 Kdet Kch Rlambda hl0 Rrec alphabar
global pabs R area volume QtorrLit Qmolec Kpump scat_Xsec
ee=1.6022E-19; %electron charge
me=9.1095E-31; % mass of electron (kg)
MO=2.6787702e-26; % mass of an Oxygen atom
MO2=2*MO; % mass of an Oxygen Molecule
Tgas = 600;
Tg=(8.61738*10^-5)*Tgas; % 600K in volts
Ti=Tg; %assumes ion temperature is equal to gas temperature
gammaO2m=0.007 % wall recombination rate of meta-stable Oxygen

%Reactor Dimensions
R=0.08 % reactor radius
L=0.075 % reactor length
l_p=L/2; % half length of reactor
area=2*pi*R*(R+L); % total surface area
volume=pi*R^2*L; % reactor volume

%Collisional energy loss per electron-ion pair created
Efactor_O=2+0.5*(1+log(MO/(2*pi*me))); % (E_e+E_i_O)/Te 
Efactor_O2=2+0.5*(1+log(MO2/(2*pi*me))); % (E_e+E_i_O2)/Te

% Loading all cross-section data to calculate Kel & Ec
load o2cross.txt -ASCII;
EnergyO2 = o2cross(:,2);
sigmaO2 = o2cross(:,3);
sigO2 = sigmaO2 * 1e-20;
load Ocross.txt -ASCII;
EnergyO = Ocross(:,1);
sigmaO = Ocross(:,2);
sigO = sigmaO * 1e-20;
scat_Xsec=7.5e-19; % elastic scattering cross-section for Oxygen in m^2
% power input
Pabs=2000; % total absorbed power in watts [adjustable] 

pabs=Pabs/(ee*volume); 

% starting pressures in mTorr (180W)

ppvec=[6.02 13.06 21.81 30.60 34.86 43.64 52.09 60.75]
Qvec=[40 50 50 50 50 50 50 50] % O2 flow rate in sccm
% atomic oxygen surface recombination rate
gammaOvec=[0.5 0.43 0.33 0.27 0.23 0.2 0.15 0.13] %gets to smaller values as pressure increases
allresults=zeros(35,length(ppvec)); % number of items to save=35
% heavy particle reaction rates
Krec=2.6E-14*sqrt(0.026/Tg);
Krec2=2.6E-14*sqrt(0.026/Tg);
Krec3=4.0E-14*sqrt(0.026/Tg);
Krec4=3.3E-17;
Kdet=1.6E-16;
Kch=2.0E-17*sqrt(0.026/Tg);


global test
vdotVector = sym(determineVdot())
test = matlabFunction(vdotVector)
%starting and end time
t0=0;
tf=90;
for ii=1:1
gammaO=gammaOvec(ii)
p=ppvec(ii);
Qsccm=Qvec(ii);
QtorrLit=Qsccm/79.05; % sccm to Torr-Liter/sec
Qmolec=4.483e17*Qsccm; % sccm to molecules/sec
Kpump=2*QtorrLit/(p*volume); % Pumping Rate coefficient (pressure in mtorr)
ng0=3.3E19*p*0.026/Tg; % m^-3
ng0_cm=ng0*1e-6 % cm^-3
%initial particle densities
nO2plusbar0=2E16
nOplusbar0=1E16
nOminusbar0=3E15
nObar0=2E17
nO2mbar0=0.01*ng0
nO20=ng0-nObar0-nO2mbar0
Te0=2;
pe0=1.5*(nO2plusbar0+nOplusbar0-nOminusbar0)*Te0; %quasineutrality condition for energy balance
v0=[nO20 nO2plusbar0 nOplusbar0 nOminusbar0 nObar0 nO2mbar0 pe0];
[t,v] = ode23s(@(t,v)dndt_MJC(t,v),[t0 tf],v0);
%[t,v] = ode23s(@dndt,[t0 tf],v0);
end
compare(:,1) = myV;
compare(:,2) = vdot;




