function power = PowerPlotOneSided(w,basic,T,decision,information)

if nargin < 5
    information = 'Bayesian';
end
if nargin < 4
    decision = 'Bayesian';
end

sigmaX = basic.sigmaX;
alpha = basic.alpha;
r = basic.r;
IN = basic.IN;
IS = basic.IS;
pN = basic.pN;
P = basic.P;
n0 = basic.n0;
mu0 = basic.mu0;
CTdiscrate = basic.CTdiscrate;
zeta = basic.zeta;

if strcmp(decision,'Frequentist')
    alphaN = @(T) sigmaX./sqrt(r.*T./2) .* norminv(1 - alpha);
    alphaS = @(T) sigmaX./sqrt(r.*T./2) .* norminv(1 - alpha);
elseif strcmp(decision,'Bayesian')
    alphaN = @(T) sigmaX./sqrt(n0 + r.*T./2) .* norminv(1 - alpha);
    alphaS = @(T) sigmaX./sqrt(n0 + r.*T./2) .* norminv(1 - alpha);
elseif strcmp(decision,'Optimal')
    if CTdiscrate > 0
        Pdisc = @(T) zeta .* (1 - exp(CTdiscrate .* P(T)./zeta))/CTdiscrate;
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
else
    error('decision mis-specified. Choose Frequentist, Bayesian, or Optimal.')
end

if strcmp(information,'Frequentist')
    power = 1 - normcdf(sqrt(r.*T)./sigmaX.*(alphaN(T) - w));
elseif strcmp(information,'Bayesian')
    if T == 0
        power = (mu0 > alphaN(T)).*ones(size(w));
    else
        power = 1 - normcdf(sqrt(r.*T./2)./sigmaX.*(((n0 + r.*T./2).*alphaN(T) - n0.*mu0)./(r.*T./2) - w));
    end
else
    error('information mis-specified. Choose Frequentist or Bayesian.')
end


