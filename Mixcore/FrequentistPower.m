function power = FrequentistPower(basic,SampleSize)
% This function calculates the power for the model specified by basic if
% the number of pairwise allocations is set to SampleSize
% 
% Revision: AA 31-10-2017

N = SampleSize;
delta = basic.relevantdifference;
sigmaX = basic.sigmaX;
alpha = basic.alpha;

effect = delta./sigmaX;
power = normcdf(sqrt(N).*effect - norminv(1-alpha./2));