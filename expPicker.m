%implement sequential experimental design picker

clear all
global expConditions
%global data
global noUnknowns
global numberOfIons
global numberOfRadicals
global noExperimentalMetrics
global levelSetUnknowns
global numberOfMaterials
global priorLB
global priorUB

priorLB = [0 0 1 1 0 0 0 0 .8 0 40 1 .001];
priorUB = [2 2 8 8 1 1 1 1 5 1 80 10 .1];
noExperimentalMetrics = 2;
numberOfIons = 2;
numberOfRadicals = 2;
levelSetUnknowns = 2;
numberOfMaterials = 1;
noUnknowns = numberOfIons*2+numberOfRadicals*2+levelSetUnknowns+1;

load allSynExpConditions
numberOfExperiments = 10;
r = randi([1,length(allSynExpConditions)],numberOfExperiments,1)
expConditions = allSynExpConditions(r,:);

%make fake real parameters for novel resist material "Algernon"
A1 = [.9 2];
B1 = [3.1 5.7];
rxnProb1 = [0.4 .7];
SCl1 = [.4 .7];
lambda1 = 1;
lambda2 = .6;
epsS1 = 50;
epsD1 = 8;
noise = .02;

A2 = [.2 .1];
B2 = [1 1];
rxnProb2 = [0.4 .8];
SCl2 = [.9 .9];
epsS2 = 11;
epsD2 = 2;

material1 = [A1 B1 rxnProb1 SCl1 lambda1 lambda2 epsS1 epsD1 noise];
material2 = [A2 B2 rxnProb2 SCl2 lambda1 lambda2 epsS2 epsD2 noise];

algernonUnknowns = [material2]
%---------------------------------------------------------------------------
% create options structure
[g, options] = createOptionsDefault();

% set case for velocities
options.etchShape = 'varyAlpha';
%---------------------------------------------------------------------------
% Modify the grid.
g.dim = 2;
g.min = [-1; -1];
g.dx = 1 / 100;
g.max = [1;1];
g.bdry = @addGhostExtrapolate;

g = processGrid(g);
%---------------------------------------------------------------------------
% create geometry, layers and mask

% set to zero for vertical stacks
options.horizontal = 1;

% get layer boundaries
% boundaries must be set in getMaterialMap. Haven't figured out an easier
% way yet
[startMatInd, endMatInd, ~, options.map] = getMaterialMap(g, options);
options.colors = getColors(options.map);
options.layer_boundaries = cat(1,startMatInd,endMatInd);
%---------------------------------------------------------------------------
% Plotting parameters
% To plot or not to plot
options.doDisplay = 0;
% Delete previous plot before showing next?
options.deleteLastPlot = 0;
% Plot in separate subplots (set deleteLastPlot = 0 in this case)?
options.useSubplots = 1;

%---------------------------------------------------------------------------
% Integration parameters.
options.tMax = 1;                  % End time.
options.plotSteps = 4;               % How many intermediate plots to produce?
options.t0 = 0;                      % Start time.
options.singleStep = 0;              % Plot at each timestep (overrides tPlot).

% Period at which intermediate plots should be produced.
options.tPlot = (options.tMax - options.t0) / (options.plotSteps - 1);
%---------------------------------------------------------------------------

err = .02; %equivalent of noise in genSynData
nOut = 10;
nIn = nOut;

for draw = 1:nIn
    for k=1:length(priorLB)-1
        theta(draw,k) = (priorUB(k)-priorLB(k)).*rand(1,1) + priorLB(k);
    end
end

%evaluate all experimental conditions for each theta
for expNo = 1:numberOfExperiments
    for j= 1:nIn
        [data, g, width(expNo,j), height,sidewallAngle, vis] = etchingWithOptions(g,options,theta(j,:),expNo);
    end
end
expObjective = width;      
U = zeros(length(expConditions),1);
for expNo = 1:length(expConditions)
    expNo
    outerSum = 0;
    for i=1:nOut
        yi = expObjective(expNo,i) + normrnd(0,err);
        innerMCSum = 0;
        for j= 1:nIn
            yj = expObjective(expNo,j) + normrnd(0,err);
            innerMCSum = normpdf(yi,yj,err)+innerMCSum;
        end
        innerTerm = innerMCSum/nIn;
        outerSum = log(normpdf(yi,expObjective(expNo,i),err)+10e-20)-log(innerTerm +10e-20) + outerSum;
    end
    U(expNo) = outerSum/nOut;
end
figure;
plot(d,U);