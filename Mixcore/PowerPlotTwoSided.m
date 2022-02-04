function power = PowerPlotTwoSided(w,basic,T,r,decision,information)

if nargin < 6
    information = 'Bayesian';
end
if nargin < 5
    decision = 'Optimal';
end
if nargin < 4 || isempty(r)
    r = basic.r;
end

sigmaX = basic.sigmaX;
alpha = basic.alpha;
IN = basic.IN;
IS = basic.IS;
pN = basic.pN;
P = basic.P;
n0 = basic.n0;
mu0 = basic.mu0;
CTdiscrate = basic.CTdiscrate;
zeta = basic.zeta;

if strcmp(decision,'Frequentist')
    alphaN = @(T) sigmaX./sqrt(r.*T./2) .* norminv(1 - alpha./2);
    alphaS = @(T) sigmaX./sqrt(r.*T./2) .* norminv(1 - alpha./2);
elseif strcmp(decision,'Bayesian')
    alphaN = @(T) sigmaX./sqrt(n0 + r.*T./2) .* norminv(1 - alpha./2);
    alphaS = @(T) sigmaX./sqrt(n0 + r.*T./2) .* norminv(1 - alpha./2);
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
    power1 = 1 - normcdf(sqrt(r.*T./2)./sigmaX.*(alphaN(T) - w));
    power2 = 1 - normcdf(sqrt(r.*T/2)./sigmaX.*(alphaS(T) + w));
    power = power1 + power2;
elseif strcmp(information,'Bayesian')
    if T == 0
        power1 = (mu0 > alphaN(T)).*ones(size(w));
        power2 = (mu0 < -alphaS(T)).*ones(size(w));
        power = power1 + power2;
    else
        power1 = 1 - normcdf(sqrt(r.*T./2)./sigmaX.*(((n0 + r.*T./2).*alphaN(T) - n0.*mu0)./(r.*T./2) - w));
        power2 =  normcdf(sqrt(r.*T./2)./sigmaX.*((-(n0 + r.*T./2).*alphaS(T) - n0.*mu0)./(r.*T./2) - w));
        power = power1 + power2;
    end
else
    error('information mis-specified. Choose Frequentist or Bayesian.')
end

