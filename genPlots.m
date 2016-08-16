function F = genPlots(theta) 

global expConditions
global data
global numberOfIons
global numberOfRadicals
global expConditions
global noUnknowns
global y
global noTrainingData
global noTestData
global algernonEtchRatesWithNoise
load allSynExpConditions
% numberOfIons = 2;
% numberOfRadicals = 1;
% noUnknowns = numberOfIons*2+numberOfRadicals*2;
% 
% %Bayesian model method
% load allSynExpConditions
% load algernonEtchRates
% 
% %Make set of training/test data %no added noise
% noTrainingData = 6;
% noTestData = 10;
% k = noTrainingData + noTestData;
% y = randsample(length(algernonEtchRates),k);

% load theta
meanInferred = mean(theta);
sdInferred = std(theta);
for k=1:length(data)
        meanEtchRate(k) = plasma(meanInferred,k);
        sdaboveEtchRate(k) = plasma(meanInferred+sdInferred,k);
        sdbelowEtchRate(k) = plasma(meanInferred-sdInferred,k);
end
exp = linspace(1,length(data),length(data));
figure;
scatter(exp,meanEtchRate)
hold on
scatter(exp,data,'r');
title('Real versus predicted (based on mean parameter values');
xlabel('Exp No');
ylabel('Etch Rate (nm/min)');

figure;
scatter(exp,sdaboveEtchRate)
hold on
scatter(exp,data,'r');
title('Real versus predicted (sd +)');
xlabel('Exp No');
ylabel('Etch Rate (nm/min)');

figure;
scatter(exp,sdbelowEtchRate)
hold on
scatter(exp,data,'r');
title('Real versus predicted (sd -');
xlabel('Exp No');
ylabel('Etch Rate (nm/min)');

data = algernonEtchRatesWithNoise(y(noTrainingData+1:length(y)));
expConditions = allSynExpConditions(y(noTrainingData+1:end),:);

for k=1:length(data)
        meanEtchRate(k) = plasma(meanInferred,k);
        sdaboveEtchRate(k) = plasma(meanInferred+sdInferred,k);
        sdbelowEtchRate(k) = plasma(meanInferred-sdInferred,k);
end

figure;
scatter(exp,meanEtchRate)
hold on
scatter(exp,data,'r');
title('Real versus predicted (Test Data - mean)');
xlabel('Exp No');
ylabel('Etch Rate (nm/min)');

figure;
scatter(exp,sdaboveEtchRate)
hold on
scatter(exp,data,'r');
title('Real versus predicted (Test Data = sd+');
xlabel('Exp No');
ylabel('Etch Rate (nm/min)');

figure;
scatter(exp,sdbelowEtchRate)
hold on
scatter(exp,data,'r');
title('Real versus predicted (Test Data = sd-');
xlabel('Exp No');
ylabel('Etch Rate (nm/min)');