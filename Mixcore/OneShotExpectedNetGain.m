function ENG = OneShotExpectedNetGain(basic,T,r)
% This function evaluates one-shot expected net gain of the model specified
% by basic at T and r.
%
% Revision: AA 28/6/2021

% Take parameters from the structure basic
c = basic.c;
pN = basic.pN;
mu0 = basic.mu0;
IS = basic.IS;
IN = basic.IN;
zeta = basic.zeta;
Tmax = basic.Tmax;
rmax = basic.rmax;
Delta = basic.Delta;
n0 = basic.n0;
sigmaX = basic.sigmaX;
cT = basic.cT;
P = basic.P;
online = basic.online;
CTdiscrate = basic.CTdiscrate;
UnkVariance = basic.UnkVariance;
xi0 = basic.UnkVarianceShape;
chi0 = basic.sigmaX^2 *(xi0 - 1);
ccap = basic.ccap;


% At the moment, this function only takes scalar values of r. We could
% extend the code to output a meshgrid of points when r is a vector.
if nargin < 2
    T = basic.T;
end
if nargin < 3
    r = basic.r;
end
if r == 0
    T = zeros(size(T));
elseif r < 0 || r > rmax
    T = -ones(size(T));
end
    
    
if CTdiscrate > 0
    Teff = @(T) (1 - exp(-CTdiscrate .* T))/CTdiscrate;
    Pdisc = @(T) zeta .* (1 - exp(-CTdiscrate .* P(T)./zeta))./CTdiscrate;
else
    Teff = @(T) T;
    Pdisc = @(T) P(T);
end

if ~UnkVariance % known sampling variance
    sigmaZ = @(T) sqrt(sigmaX^2 .* r .* T ./ (n0 .* (2.*n0 + r .* T))); % standard deviation of the predictive mean
    % alphaN and alphaS are the costs of adoption of technology N and S, resp.
    % We get NaN when IN=0 and either P(T)=0 or pN=0 because we are dividing 0
    % by 0. This chunk of code removes that issue by setting alphaN =-mu0 when
    % IN=0 and setting alphaS=mu0 when IS=0.
    if IN ==0
        alphaN = @(T) 0;
    else
        alphaN = @(T) IN ./ ((1-pN) .* Pdisc(T));
    end
    if IS == 0
        alphaS = @(T) 0;
    else
        alphaS = @(T) IS ./ (pN .* Pdisc(T));
    end
    
    % The values of T are divided in three regions. T=0, 0<T<H-\Delta, T out of
    % bounds. Tzeros are the logical indices where T=0, Toob are the indices
    % where T is out of bounds, and Tinterior are the logical indices where
    % 0<T<H-\Delta
    Tzeros = T == 0;
    Toob = T < 0 | T > Tmax;
    Tinterior = T > 0 & T <= Tmax;
    % Define ENG for the three regions
    ENG = zeros(size(T));
    ENG(Tzeros) = max([0 , Pdisc(0).*(1-pN).*mu0 - IN , -Pdisc(0).*pN.*mu0 - IS]);
    ENG(Toob) = -inf;
    T = T(Tinterior);
    ENG(Tinterior) = 0.5 .* online .* r .* Teff(T) .* (1 - 2.*pN) .* mu0 ...
        + exp(-CTdiscrate.*(T + Delta)).*Pdisc(T).*(1 - pN).*sigmaZ(T).*Psi((alphaN(T) - mu0)./sigmaZ(T))...
        + exp(-CTdiscrate.*(T + Delta)).*Pdisc(T).*pN.*sigmaZ(T).*Psi((alphaS(T) + mu0)./sigmaZ(T))...
        - ccap(r) - c.*r.* Teff(T) - cT.*Teff(T);
    
else % If the variance is unknown
    sigmaZ = @(T) sqrt(chi0 .* r .* T ./ (xi0 .* n0 .* (2.*n0 + r .* T))); % standard deviation of the predictive mean
    % alphaN and alphaS are the costs of adoption of technology N and S, resp.
    % We get NaN when IN=0 and either P(T)=0 or pN=0 because we are dividing 0
    % by 0. This chunk of code removes that issue by setting alphaN =-mu0 when
    % IN=0 and setting alphaS=mu0 when IS=0.
    if IN ==0
        alphaN = @(T) 0;
    else
        alphaN = @(T) IN ./ ((1-pN) .* Pdisc(T));
    end
    if IS == 0
        alphaS = @(T) 0;
    else
        alphaS = @(T) IS ./ (pN .* Pdisc(T));
    end
    
    % The values of T are divided in three regions. T=0, 0<T<H-\Delta, T out of
    % bounds. Tzeros are the logical indices where T=0, Toob are the indices
    % where T is out of bounds, and Tinterior are the logical indices where
    % 0<T<H-\Delta
    Tzeros = T == 0;
    Toob = T < 0 | T > Tmax;
    Tinterior = T > 0 & T <= Tmax;
    % Define ENG for the three regions
    ENG = zeros(size(T));
    ENG(Tzeros) = max([0 , Pdisc(0).*(1-pN).*mu0 - IN , -Pdisc(0).*pN*mu0 - IS]);
    ENG(Toob) = -inf;
    T = T(Tinterior);
    ENG(Tinterior) = 0.5 .* online .* r .* Teff(T) .* (1 - 2.*pN) .* mu0 ...
        + exp(-CTdiscrate.*(T + Delta)).*Pdisc(T).*(1 - pN).*sigmaZ(T).*tPsi((alphaN(T) - mu0)./sigmaZ(T), 2.*xi0)...
        + exp(-CTdiscrate.*(T + Delta)).*Pdisc(T).*pN.*sigmaZ(T).*tPsi((alphaS(T) + mu0)./sigmaZ(T), 2.*xi0)...
        - ccap(r) - c.*r.* Teff(T) - cT.*Teff(T);
    
    
end
end


function y = Psi(x)
index = x == Inf;
y = zeros(size(x));
y(index) = 0;
x = x(~index);
y(~index) = normpdf(x) - x .* (1 - normcdf(x));
end

function y = tPsi(x,nu)
index = x == Inf;
y = zeros(size(x));
y(index) = 0;
x = x(~index);
y(~index) = (nu - x.^2)./(nu - 1) .* tpdf(x,nu) - x.*(1 - tcdf(x,nu));
end