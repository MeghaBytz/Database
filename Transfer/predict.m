clear all
close all

%this method uses Bayesian inference to infer ion sputtering yield,
%reaction probabilty of neutrals and ion stimulated desorption (all
%hard to measure constants in the surface kinetics model for etching). 

%declare global variables
global proposalLB
global proposalUB
global data
global noUnknowns
global priorLB
global priorUB
global expRows
global expColumns
global noExperiments
global proposalMean
global proposalSD
global numberOfIons
global numberOfRadicals
global expConditions
global y
global noTrainingData
global noTestData
global noExperimentalMetrics
%global algernonEtchRatesWithNoise

noExperimentalMetrics = 2;
numberOfIons = 2;
numberOfRadicals = 1;
noUnknowns = numberOfIons*2+numberOfRadicals*2+1;
proposalLB = [0 0 0 0 0 0 .01]
proposalUB = [1 1 1 1 1 1 1];
proposalMean = [5 5 5 5 5 5 7.5];
proposalSD =  ones(noUnknowns,1)*0.5;
priorLB = [0 0 0 0 0 0 .01]
priorUB = [10 10 10 10 10 10 15];


%Bayesian model method
load allSynExpConditions
load algernonEtchRates

%add noise to data
sd = 7;
%algernonEtchRatesWithNoise = algernonEtchRates + sd*randn(size(algernonEtchRates));
algernonEtchRatesWithNoise = [22.2165602702964,77.6055590499576,236.575582379589,450.079040977680,35.3605596462460,70.6079924218272,245.367890373969,447.545239908143,25.0484496346538,70.6440344816242,238.695051213768,439.411022589560,15.7321038990669,31.7014270882974,115.828400282430,237.469076616231,10.5454450551198,53.6399633481009,120.061838020585,236.961835578835,15.5285213002027,40.0578309474414,128.020131273519,236.206717710504,15.9208596916653,17.1430241490525,81.1144309189445,157.402620742333,14.4390304437448,20.8172399582352,88.3083270967013,133.897334316689,9.82944786941629,22.1906680675970,80.0875220484532,139.560246665357,18.8823335115934,13.5685309353505,60.4025995287640,110.386880988572,1.47110312941420,22.2722872825582,55.8669645567098,93.3237387216162,1.00731675432532,16.0682864511765,60.6924137091444,103.849318871479,7.65333337435853,15.1181668359758,47.4174480095480,82.8128201515366,7.74915868620683,22.5047721584586,51.5817305991003,83.5107479164520,18.3985706506380,9.43517827473717,57.0696300898955,78.1852597843657,6.96798216633548,12.3161199359238,33.4681493441422,66.7214465758618,16.3399185037639,9.73238562633700,37.7724261665344,73.6603826623381,0.728989378135380,10.5091615946335,40.9625756896852,61.2621312774054,9.91938808038657,12.8939391004239,42.5867233604335,41.1948034140787,7.15680329469013,6.39041390401875,31.3406350570565,58.9411969091671,5.63879861664421,10.4376739081442,30.8180685613177,54.1886221470530,19.5671017424484,6.84846336034813,39.1141920572799,36.7292380222424,16.0342656498202,0.579706280545058,29.9336121581765,49.6670931327384,-5.31825697617586,4.88054899988278,23.2380803703649,44.4425763342825]

%Make set of training/test data %no added noise
noTrainingData = 3;
noTestData = 10;
k = noTrainingData + noTestData;
y = randsample(length(algernonEtchRatesWithNoise),k);
%data = algernonEtchRatesWithNoise(y(1:noTrainingData));
expConditions = allSynExpConditions(y(1:noTrainingData),:);
data = [.5 1; .2 .9; .3 .91]
%realParameters = [.05 .05 .05 .05 0.25]
%Make initial guess for unknown parameters
current = ones(noUnknowns,1);
Likelihood(current)
%Peform MH iterations

N = 10;
theta = zeros(N*noUnknowns,noUnknowns);
acc = zeros(1,noUnknowns);
PosteriorCurrent = Posterior(current,1);
% for i = 1:burnin    % First make the burn-in stage
%     for j=1:noUnknowns
%       [alpha,t, a,prob, PosteriorCatch] = MetropolisHastings(theta,current,PosteriorCurrent,j);
%     end
% end
for cycle = 1:N  % Cycle to the number of samples
    ind = 1;
     for m=1:noUnknowns % Cycle to make the thinning
            [alpha,t, a,prob, PosteriorCatch] = MetropolisHastings(theta,current,PosteriorCurrent,m);
            theta((cycle-1)*noUnknowns+m,:) = t;        % Samples accepted
            AlphaSet(cycle,m) = alpha;
            current = t;
            PosteriorCurrent = PosteriorCatch;
            acc(m) = acc(m) + a;  % Accepted ?
     end
end
 

accrate = acc/N;     % Acceptance rate,. 

%calculate final etch rate from Bayes using mean values
for i =1:noUnknowns
    inferred(i) = mean(theta(:,i));
end
for k=1:length(data)
        inferredEtchRate(k) = plasma(inferred,k);
end

%Compare real versus prediction
exp = linspace(1,length(data),length(data));
scatter(exp,inferredEtchRate)
hold on
scatter(exp,data,'r');
title('Real versus predicted (based on mean parameter values');
xlabel('Exp No');
ylabel('Etch Rate (nm/min)');
figure;
    for i =1:noUnknowns
        subplot(3,3,i);
        outputTitle = sprintf('Unknown %d',i);
        hist(theta(:,i)); 
 %       f.Normalization = 'probability';
        %counts = f.Values;
%         hold on
%         line([real(i) real(i)],[0 max(counts)], 'Color', 'r')
%         hold on
%         xmin = min(theta(:,i));
%         xmax = max(theta(:,i));
%         x = xmin:1:xmax;
%         hold on
%         line([priorLB(i) priorLB(i)],[0 max(counts)], 'Color', 'g')
%         hold on
%         line([priorUB(i) priorUB(i)],[0 max(counts)], 'Color', 'g')
        title(outputTitle);
        xlabel('Value');
        ylabel('Frequency');
    end



figure; 
for k =1:noUnknowns
    outputTitle = sprintf('Unknown %d',k);
    subplot(3,3,k);
    plot(theta(:,k));
    hold on
    line([0 N],[real(k) real(k)], 'Color', 'r')
    hold on
    line([0 N],[proposalLB(i) proposalLB(i)], 'Color', 'g')
    line([0 N],[proposalUB(i) proposalUB(i)], 'Color', 'g')
    title(outputTitle);
    xlabel('Cycle #');
    ylabel('Value');
end

%genPlots(theta)
