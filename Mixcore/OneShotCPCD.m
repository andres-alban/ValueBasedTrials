function CPCD = OneShotCPCD(w,basic,T,r)
% This function computes the conditional probability of correct selection
% (CPCS) for a trial defined by basic with decision variables T and r, for
% the true INMB given by vector w.

sigmaX = basic.sigmaX;
alpha = basic.alpha;
IN = basic.IN;
IS = basic.IS;
pN = basic.pN;
P = basic.P;
n0 = basic.n0;
mu0 = basic.mu0;
delta = basic.relevantdifference;
Delta = basic.Delta;
CTdiscrate = basic.CTdiscrate;
zeta = basic.zeta;

if nargin<4
    r = basic.r;
end

if CTdiscrate > 0
    Pdisc = @(T) zeta .* (1 - exp(-CTdiscrate .* P(T)./zeta))./CTdiscrate;
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
    alphaS = @(T) IS ./ (pN .* Pdisc(T));
end


IndexS = w < -alphaS(T);
IndexN = w > alphaN(T);
IndexM = ~(IndexS | IndexN);

CPCD = zeros(size(w));
% The following ignores that w=alphaN or alphaS requires special treatment
% because any choice is correct.
if r*T <= 0
    CPCD(IndexS) = (mu0 < -alphaS(T));
    CPCD(IndexN) = (mu0 > alphaN(T));
    CPCD(IndexM) = (mu0 <= alphaN(T)) & (mu0 >= -alphaS(T));
else
    CPCD(IndexS) = 1 - normcdf((2.*n0.*(alphaS(T) + mu0) + r.*T.*(alphaS(T) + w(IndexS)))./(sigmaX.*sqrt(2.*r.*T)));
    CPCD(IndexN) = 1 - normcdf((2.*n0.*(alphaN(T) - mu0) + r.*T.*(alphaN(T) - w(IndexN)))./(sigmaX.*sqrt(2.*r.*T)));
    CPCD(IndexM) = normcdf((2.*n0.*(alphaS(T) + mu0) + r.*T.*(alphaS(T) + w(IndexM)))./(sigmaX.*sqrt(2.*r.*T)))...
        + normcdf((2.*n0.*(alphaN(T) - mu0) + r.*T.*(alphaN(T) - w(IndexM)))./(sigmaX.*sqrt(2.*r.*T))) - 1;
end




