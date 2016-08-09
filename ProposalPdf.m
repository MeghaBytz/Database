function F = ProposalPdf(current,new,index,chainSD)
proposalPDFTime = tic;
%use current instead of center to change back to random walk
global proposalLambda;
global proposalSD
global proposalRecord
global proposalMean
global proposalLB
global proposalUB
q = 1;
q=1;
q = lognpdf(new(index),proposalMean(index),chainSD)*q;
F = log(q);
end