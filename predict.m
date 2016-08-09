clear all
close all

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


numberOfIons = 2;
numberOfRadicals = 1;
noUnknowns = numberOfIons*2+numberOfRadicals;
proposalLB = zeros(noUnknowns,1);
proposalUB = [10 10 10 10 1];
proposalMean = [5 5 5 5 .5]
proposalSD =  [2 2 2 2 1];
priorLB = zeros(noUnknowns,1);
priorUB = [10 10 10 10 1];
data = [3 5 3];

%Bayesian model method
load expConditions
%Make initial guess for unknown parameters
current = zeros(noUnknowns,1);
%Peform MH iterations
PosteriorCurrent = Posterior(current,1);
N = 100;
theta = zeros(N*noUnknowns,noUnknowns);
acc = zeros(1,noUnknowns);
for i = 1:noUnknowns
          current = ProposalFunction(theta,current,i);
end
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

figure;
    for i =1:noUnknowns
        subplot(3,3,i);
        outputTitle = sprintf('Unknown %d',i);
        paramEdges = [0:2:proposalUB(i)];
        hist(theta(:,i),paramEdges);
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

