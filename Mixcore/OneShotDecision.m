function [DN,DS,DM] = OneShotDecision(basic,T)
% This function computes the probability of adopting each alternative for a
% trial defined by basic with recruitment rate basic.r and duration T.

% Take parameters from the structure basic
c = basic.c;
pN = basic.pN;
mu0 = basic.mu0;
IS = basic.IS;
IN = basic.IN;
Delta = basic.Delta;
n0 = basic.n0;
sigmaX = basic.sigmaX;
cfix = basic.cfix;
csite = basic.csite;
s = basic.s;
r = basic.r;
P = basic.P;
online = basic.online;
CTdiscrate = basic.CTdiscrate;


% Define useful variables that are also defined in the paper
sigmaZ = @(T) sqrt(sigmaX^2 .* r .* T ./ (n0 .* (n0 + r .* T))); % standard deviation of the predictive mean

if CTdiscrate > 0
    Pdisc = @(T) rho .* (1 - exp(CTdiscrate .* P(T)./rho)/CTdiscrate);
else
    Pdisc = @(T) P(T);
end
% alphaN and alphaS are the financial per person costs of adoption of
% technology N and S, resp. We get NaN when IN=0 and either P(T)=0 or
% pN=0 because we are dividing 0 by 0. This chunk of code removes that
% issue.
if IN == 0
    alphaN = @(T) 0;
else
    alphaN = @(T) IN ./ ((1-pN) .* Pdisc(T));
end
if IS == 0
    alphaS = @(T) 0;
else
    alphaS = @(T) -IS ./ (pN .* Pdisc(T));
end

DN = 1 - normcdf((alphaN(T) - mu0)./sigmaZ(T));
DS = 1 - normcdf((alphaS(T) + mu0)./sigmaZ(T));
DM = 1 - DN - DS;



