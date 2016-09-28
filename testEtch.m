clear all
close all

global data
global expConditions


load y
load theta_withNoise
load('algernonEtchRatesWithNoise.mat')
test = algernonEtchRatesWithNoise
% testSet = y;
% data = algernonEtchRatesWithNoise(testSet(7:end));
% expConditions = allSynExpConditions(testSet(7:end),:);
% 
% inferred = mean(theta_withNoise);
% 
% for k=1:length(data)
%         inferredTestEtchRates(k) = plasma(inferred,k);
% end