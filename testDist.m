% maxdist - Compute probabilities from Maxwell distribution
clear all;        % Clear memory
close all
%@ Initialize variables
m = 4.48e-26;     % Mass of nitrogen molecule
k = 1.38e-23;     % Boltzmann's constant
ee = 1.602*10^-19;
Te = .052; %(volts)
EvtoK = 1/(8.621738e-5);
Te = Te*EvtoK;
draws = 1000;
drawDistribution1 = [];
for k = 1:draws
    x1 = rand;
    x2 = rand;
    y1 = sqrt(k*Te/(m))*sqrt(-2*log(x1))*sin(2*pi*x2);
    y2 = sqrt(k*Te/(m))*sqrt(-2*log(x1))*sin(2*pi*x2);
    drawDistribution1 = [drawDistribution1 y1 y2]
end
figure(1)
hist(drawDistribution1)
std(drawDistribution1)
vi = 
a = sqrt(k*Te/(2*m))