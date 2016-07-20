myV(1)=-k1 * n1 * n8 + -k2 * n1 * n8 + -k3 * n1 * n8 + -k4 * n1 * n8 + -k5 * n1 * n8 + -k6 * n1 * n8 + k9 * n2 * n5 + -k10 * n1 * n4 + -k11 * n1 * n8 + k12 * n6 * n8 + k16 * n5 * n6 - Kpump*n1 + k17 * np3 * np5 * Vrec/V + k20/2 * np2 * np9 + k21 * np6 * np9 + k22 * np3 * np9 + Q1/V;
myV(2)=k2 * n1 * n8 + k3 * n1 * n8 + 2*k5 * n1 * n8 + k6 * n1 * n8 + -k7 * n2 * n8 + 2*k8 * n3 * n8 + -k9 * n2 * n5 + k10 * n1 * n4 + k14 * n6 * n8 + 2*k15 * n6 * n8 + k16 * n5 * n6 - Kpump*n2 + k17 * np3 * np5 * Vrec/V + 3*k18 * np3 * np5 * Vrec/V + 2*k19 * np4 * np5 * Vrec/V + -k20 * np2 * np9 + k23 * np4 * np9;
myV(3)=k1 * n1 * n8 + -k8 * n3 * n8 + k10 * n1 * n4 + k13 * n6 * n8 + -k17 * np3 * np5 * Vrec/V + -k18 * np3 * np5 * Vrec/V + -k22 * np3 * np9;
myV(4)=k3 * n1 * n8 + k4 * n1 * n8 + k7 * n2 * n8 + -k10 * n1 * n4 + -k19 * np4 * np5 * Vrec/V + -k23 * np4 * np9;
myV(5)=k2 * n1 * n8 + k4 * n1 * n8 + -k9 * n2 * n5 + k14 * n6 * n8 + -k16 * n5 * n6 + -k17 * np3 * np5 * Vrec/V + -k18 * np3 * np5 * Vrec/V + -k19 * np4 * np5 * Vrec/V;
myV(6)=k11 * n1 * n8 + -k12 * n6 * n8 + -k13 * n6 * n8 + -k14 * n6 * n8 + -k15 * n6 * n8 + -k16 * n5 * n6 - Kpump*n6 + -k21 * np6 * np9;
myV(7)=k6 * n1 * n8 - Kpump*n7;
myV(8)=pabs - Ec_1*k1 * n1 * n8 - Ec_2*k7 * n2 * n8 - Ew_1*k22 * np3 * np9 - Ew_2*k23 * np4 * np9;
