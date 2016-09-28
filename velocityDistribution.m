% maxdist - Compute probabilities from Maxwell distribution
clear all;        % Clear memory
help maxdist;     % Print header
v = 30;
%@ Initialize variables
m = 4.48e-26;     % Mass of nitrogen molecule
k = 1.38e-23;   
ee = 1.602*10^-19; 
T = .052*ee;
% Compute two frequently used constants
Constant1 = 4*pi*(m/(2*pi*k*T))^(3/2);
Constant2 = m/(2*k*T);
ProbDensity = Constant1* v^2 *exp(-Constant2*v^2);

% Enter the lower and upper bounds for integration
vLow  = input('Enter lower bound speed (m/s): ');
vHigh = input('Enter upper bound speed (m/s): ');

%@ Perform running sum to estimate the integral
NSum = 300;         % Number of elements in the sum
DeltaV = (vHigh - vLow)/(NSum - 1);    % Delta v 
Probability = 0.0;  % Set running sum to zero
for iSum=1:NSum     % Loop until hitting the break
  %@ Compute probability density for speed v
  v = vLow + (iSum-1)*DeltaV;
  ProbDensity = Constant1* v^2 *exp(-Constant2*v^2);
  %@ Accumulate running sum
  Probability = Probability + ProbDensity*DeltaV;
end

%@ Set up the plotting variables
vPlotMax = 5*sqrt(k*T/m);     % Maximum speed to plot
NPoints = 100;                % Number of points used in plot
for iPoint=1:NPoints
  % Compute x-axis values for the plot
  vPlot(iPoint) = (iPoint-1)/(NPoints) * vPlotMax; 
  % Compute y-axis values for the plot
  fPlot(iPoint) = Constant1 * vPlot(iPoint)^2 *   ...
                     exp(-Constant2 * vPlot(iPoint)^2);
end
					 
%@ Create the plot
clf;  figure(gcf);       % Clear the figure; bring window forward
plot(vPlot,fPlot,'-');   % Plot the probability density function
hold on;                 % Hold the graph to add vertical dashed lines
% Compute values of probability at vLow and vHigh
fLow = Constant1 * vLow^2 * exp(-Constant2 * vLow^2);
fHigh = Constant1 * vHigh^2 * exp(-Constant2*vHigh^2);
% Draw vertical dashed lines between the lower and upper 
% bounds of integration
for v=vLow:(vHigh-vLow)/30:vHigh
 f = Constant1 * v^2 * exp(-Constant2*v^2);
 plot([v v],[0 f],'--r');
end
plot(vLow,0,'*b',vHigh,0,'*b');   % Mark endpoints on the x-axis
xlabel('Speed (m/s)');
ylabel('Distribution function');
title( sprintf('Probability in dashed area = %g',Probability) );