function vdot = governingEquations()
vdot(1)=Qmolec/volume + Krec*nO2plus*nOminus*vol_rec/volume...
 + Kdet*nOminusbar*nObar + Kdeex*nO2mbar*ne0...
 + Krec4*nO2mbar*nOminusbar + Kion*nO2plus + KO2m*nO2m...
 + 0.5*KO*nO - (Kiz1 + Katt + Kdiss + Kiz3 + Kiz4 + Kex)*nO2*ne0...
 - Kch*nOplusbar*nO2 - Kpump*nO2;
% for nO2plus
vdot(2)=Kiz1*nO2*ne0 + Kizm*nO2mbar*ne0 + Kch*nOplusbar*nO2...
 -(Krec + Krec2)*nO2plus*nOminus*vol_rec/volume...
 - Kei*ne0*nO2plusbar - Kion*nO2plus;
% for nOplus
vdot(3)=Kiz2*nObar*ne0 + (Kiz3 + Kiz4)*nO2*ne0...
 - Krec3*nOplus*nOminus*vol_rec/volume - Kch*nOplusbar*nO2...
 - Kion*nOplus;
% for nOminus
vdot(4)=(Katt+Kiz3)*nO2*ne0 + Kattm*nO2mbar*ne0...
 - ((Krec+Krec2)*nO2plus*nOminus+Krec3*nOplus*nOminus)...
 *vol_rec/volume - Kdet*nOminusbar*nObar - Krec4*nOminusbar*nO2mbar;
% for nO
vdot(5)=2*Kei*ne0*nO2plusbar + (2*Kdiss+Katt+Kiz4)*ne0*nO2...
 + (Krec+3*Krec2)*nO2plus*nOminus*vol_rec/volume...
 + 2*Krec3*nOplus*nOminus*vol_rec/volume + Kch*nOplusbar*nO2...
 + (Kattm+2*Kdism)*nO2mbar*ne0 + Krec4*nOminusbar*nO2mbar...
 + Kion*nOplus - Kiz2*nObar*ne0 - Kdet*nOminusbar*nObar...
 - KO*nO - Kpump*nObar;
% for nO2m
vdot(6)=Kex*nO2*ne0 - (Kizm+Kattm+Kdeex+Kdism)*nO2mbar*ne0...
 - Krec4*nO2mbar*nOminusbar - KO2m*nO2m - Kpump*nO2mbar; 
% for power balance
vdot(7)=pabs - Ec_O2*Kiz1*nO2*ne0 - Ec_O*Kiz2*nObar*ne0...
 - Eei_O*Kion*nOplus - Eei_O2*Kion*nO2plus; 
end
