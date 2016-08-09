function [F, chainSD] = ProposalFunction(theta,current,index)
global proposedParameterRecord
global proposalLB
global proposalUB
global proposalMean
global proposalSD
parameter = current;
s = 2.4;
epsilon = 10e-3;
variance = s*(var(theta(:,index)))+s*epsilon;
chainSD = sqrt(log(variance/((proposalMean(index))^2)+1));
parameter(index) = lognrnd(proposalMean(index),chainSD);
F = parameter;
end
