clear all
global expConditions
global data
global noUnknowns
global numberOfIons
global numberOfRadicals

numberOfIons = 2;
numberOfRadicals = 1;
noUnknowns = numberOfIons*2+numberOfRadicals*2+1;

load allSynExpConditions
load allSynData
%Make set of training/test data
noTrainingData = 5;
noTestData = 5;
k = noTrainingData + noTestData;
y = randsample(length(allSynData),k);
data = allSynData(y(1:noTrainingData));
expConditions = allSynExpConditions;

%make fake real parameters for novel material "Algernon"
A = [0.05 .1];
B = [0.08 .3];
rxnProb = 0.7;
SCl = .8;
noise = 1;
algernonUnknowns = [A B rxnProb SCl noise];
for i = 1:length(expConditions)
 algernonEtchRates(i,:) = plasma(algernonUnknowns,i);
end