function N = FrequentistSampleSize(basic)
% This function calculates the minimum sample size for the specified type I
% error basic.alpha, power basic.power and smallest relevant difference
% basic.relevantdifference. The output is the number of pairwise
% allocations

delta = basic.relevantdifference;
sigmaX = basic.sigmaX;
alpha = basic.alpha;
power = basic.power;

N = norminv(1 - alpha/2) + norminv(power);
effect = delta./sigmaX;
N = N^2 ./ effect^2;
N = ceil(N);