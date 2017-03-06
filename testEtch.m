clear all
close all

global data
global expConditions


load y
load('theta_withNoise.mat')
load('algernonEtchRatesWithNoise.mat')
test = algernonEtchRates;
load allSynExpConditions
load algernonEtchRates
testSet = y;
data = algernonEtchRates(testSet(7:end));
expConditions = [6.02,4,0.5,500]

inferred = mean(theta);
Te = 2;
etchRates = plasma(inferred,1);