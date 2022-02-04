function ENG = OneShotExpectedNetGainFreq(basic,r,Q)
% This function evaluates one-shot expected net gain of the model specified
% by basic when the adoption decision is two-sided frequentist hypothesis 
% testing with type 1 error specified by basic.alpha
%

% Take parameters from the structure basic
c = basic.c;
pN = basic.pN;
mu0 = basic.mu0;
IS = basic.IS;
IN = basic.IN;
zeta = basic.zeta;
H = basic.H;
Tmax = basic.Tmax;
Delta = basic.Delta;
n0 = basic.n0;
sigmaX = basic.sigmaX;
cfix = basic.cfix;
cr = basic.cr;
cT = basic.cT;
P = basic.P;
online = basic.online;
CTdiscrate = basic.CTdiscrate;
UnkVariance = basic.UnkVariance;
xi0 = basic.UnkVarianceShape;
chi0 = basic.sigmaX^2 *(xi0 - 1);
alpha = basic.alpha;
ccap = basic.ccap;

if nargin < 3
    Q = FrequentistSampleSize(basic);
end
    
    
if CTdiscrate > 0
    Teff = @(T) (1 - exp(-CTdiscrate .* T))/CTdiscrate;
    Pdisc = @(T) zeta .* (1 - exp(-CTdiscrate .* P(T)./zeta))./CTdiscrate;
else
    Teff = @(T) T;
    Pdisc = @(T) P(T);
end

if ~UnkVariance % known sampling variance
    sigmaZ = @(Q) sqrt(sigmaX^2 .* Q ./ (n0 .* (2.*n0 + Q))); % standard deviation of the predictive mean
    % alphaN and alphaS are the costs of adoption of technology N and S, resp.
    % We get NaN when IN=0 and either P(T)=0 or pN=0 because we are dividing 0
    % by 0. This chunk of code removes that issue by setting alphaN =-mu0 when
    % IN=0 and setting alphaS=mu0 when IS=0.
    
    z_alphahalf = norminv(1 - alpha./2);
    alphaN = @(T) (n0.*mu0 + sqrt(Q).*sigmaX.*z_alphahalf)./(n0 + Q);
    alphaS = @(T) (-n0.*mu0 + sqrt(Q).*sigmaX.*z_alphahalf)./(n0 + Q);
    
    
    % The values of T are divided in three regions. T=0, 0<T<H-\Delta, T out of
    % bounds. Tzeros are the logical indices where T=0, Toob are the indices
    % where T is out of bounds, and Tinterior are the logical indices where
    % 0<T<H-\Delta
    rzeros = r == 0;
    roob = r < 0 | r > zeta;
    rinterior = r > 0 & r <= zeta;
    % Define ENG for the three regions
    ENG = zeros(size(r));
    ENG(rzeros) = 0;
    ENG(roob) = -inf;
    r = r(rinterior);
    ENG(rinterior) = 0.5 .* online .* r .* Teff(2.*Q./r) .* (1 - 2.*pN) .* mu0 ...
        + exp(-CTdiscrate.*(2.*Q./r + Delta)).*Pdisc(2.*Q./r).*(1 - pN).*sigmaZ(Q).*normpdf((alphaN(2.*Q./r) - mu0)./sigmaZ(Q))...
        + exp(-CTdiscrate.*(2.*Q./r + Delta)).*(Pdisc(2.*Q./r).*(1 - pN).*mu0 - IN).*(1 - normcdf((alphaN(2.*Q./r) - mu0)./sigmaZ(Q)))...
        + exp(-CTdiscrate.*(2.*Q./r + Delta)).*Pdisc(2.*Q./r).*pN.*sigmaZ(Q).*normpdf((alphaS(2.*Q./r) + mu0)./sigmaZ(Q))...
        + exp(-CTdiscrate.*(2.*Q./r + Delta)).*(Pdisc(2.*Q./r).*pN.*mu0 - IS).*(1 - normcdf((alphaS(2.*Q./r) + mu0)./sigmaZ(Q)))...
        - ccap(r) - c.*r.* Teff(2.*Q./r) - cT.*Teff(2.*Q./r);
    
else % If the variance is unknown
    error('Not available for unknown variance yet')
%     sigmaZ = @(T) sqrt(chi0 .* r .* T ./ (xi0 .* n0 .* (2.*n0 + r .* T))); % standard deviation of the predictive mean
%     % alphaN and alphaS are the costs of adoption of technology N and S, resp.
%     % We get NaN when IN=0 and either P(T)=0 or pN=0 because we are dividing 0
%     % by 0. This chunk of code removes that issue by setting alphaN =-mu0 when
%     % IN=0 and setting alphaS=mu0 when IS=0.
%     
%     z_alphahalf = tinv(1 - alpha./2,Q);
%     alphaN = @(T) (n0.*mu0 + sqrt(Q).*sigmaX.*z_alphahalf)./(n0 + Q);
%     alphaS = @(T) (-n0.*mu0 + sqrt(Q).*sigmaX.*z_alphahalf)./(n0 + Q);
%     
%     % The values of T are divided in three regions. T=0, 0<T<H-\Delta, T out of
%     % bounds. Tzeros are the logical indices where T=0, Toob are the indices
%     % where T is out of bounds, and Tinterior are the logical indices where
%     % 0<T<H-\Delta
%     rzeros = r == 0;
%     roob = r < 0 | r > zeta;
%     rinterior = r > 0 & r <= zeta;
%     % Define ENG for the three regions
%     ENG = zeros(size(r));
%     ENG(rzeros) = max([0 , P(0).*(1-pN).*mu0 - IN , -P(0).*pN*mu0 - IS]);
%     ENG(roob) = -inf;
%     r = r(rinterior);
%     ENG(rinterior) = 0.5 .* online .* r .* Teff(2.*Q./r) .* (1 - 2.*pN) .* mu0 ...
%         + exp(-CTdiscrate.*(2.*Q./r + Delta)).*Pdisc(2.*Q./r).*(1 - pN).*sigmaZ(2.*Q./r).*tPsi((alphaN(2.*Q./r) - mu0)./sigmaZ(2.*Q./r), 2.*xi0)...
%         + exp(-CTdiscrate.*(2.*Q./r + Delta)).*Pdisc(2.*Q./r).*pN.*sigmaZ(2.*Q./r).*tPsi((alphaS(2.*Q./r) + mu0)./sigmaZ(2.*Q./r), 2.*xi0)...
%         - cfix - c.*r.* Teff(2.*Q./r) - cr.*r - cT.*Teff(2.*Q./r);
%     
    
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