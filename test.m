clear 
close all

global data
global expConditions
global noUnknowns
global numberOfIons
global numberOfRadicals
global viewFactor

%Bayesian model method
load allSynExpConditions
load algernonEtchRates
numberOfIons = 2;
numberOfRadicals = 1;
noUnknowns = numberOfIons*2+numberOfRadicals*2+1;
data = algernonEtchRates;
expConditions = allSynExpConditions;
sigma = [0.588269814728046,0.523948720678330,0.110347557790930,0.0914668271114918,0.457056073946934,0.0817223074092960,1.80538265845552];
testParam = [1.02538716198995,1.07686376359793,0.128144181651286,0.576087305315432,0.685995636526662,1.36435668956040,4.71979138765046];
testParamPlusSigma = testParam + sigma;
testParamMinusSigma = testParam - sigma;
%viewFactor = radiosity();

%need to make 2d matrix
for i = 1:1%length(data)
 predEtchRates3training = plasma(testParam,i)*viewFactor(:,1);
end

% for i = 1:length(data)
%  predEtchRates3trainingPlusSigma(i) = plasma(testParamPlusSigma,i);
% end
% 
% for i = 1:length(data)
%  predEtchRates3trainingMinusSigma(i) = plasma(testParamMinusSigma,i);
% end