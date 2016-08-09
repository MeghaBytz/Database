function [F] = Likelihood(current)
global data;
global likelihoodRecord
global etchRecord
likeData = 0;
noise = 1; %change this back to unknown noise parameter eventually
for i=1:length(data)
        etchRate = plasma(current,i);
        Like = normpdf(data(i),etchRate,noise) + 10e-20;
        likeData = log(Like) + likeData;
end
likelihoodRecord = [likelihoodRecord likeData];
F = likeData;
