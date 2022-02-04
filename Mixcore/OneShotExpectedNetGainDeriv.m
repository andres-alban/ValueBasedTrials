function ENGderivative = OneShotExpectedNetGainDeriv(basic,T)
% This function is outdated. Do not use.
% This function calculates the derivative of the expected net gain with
% respect to T defined by basic evaluated at T. This is currently not
% implemented in any other function but can potentially be used for more
% efficient optimization.
%
% Revision: AA 1-11-2017

% Take parameters from basic
c = basic.c;
pN = basic.pN;
mu0 = basic.mu0;
IS = basic.IS;
IN = basic.IN;
H = basic.H;
Delta = basic.Delta;
n0 = basic.n0;
sigmaX = basic.sigmaX;
r = basic.r;
P = basic.P;
P_T = basic.P_T;
online = basic.online;
CTdiscrate = basic.CTdiscrate;
UnkVariance = basic.UnkVariance;

% Define useful variables defined in the paper
sigmaZ = @(T) sqrt(sigmaX^2 .* r .* T ./ (n0 .* (2.*n0 + r .* T)));
sigmaZ_T = @(T) 1/2 .* sqrt(n0.*sigmaX^2.*r ./ ((2.*n0 + r.*T).^3 .* T)); % Derivative of sigmaZ with respect to T

if UnkVariance
    error('This function does not evaluate the derivative for unknown variance.')
end

if CTdiscrate > 0
    Teff = @(T) (1 - exp(-CTdiscrate .* T))/CTdiscrate;
    Pdisc = @(T) rho .* (1 - exp(CTdiscrate .* P(T)./rho)/CTdiscrate);
else
    Teff = @(T) T;
    Pdisc = @(T) P(T);
end

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
    Toob = T < 0 | T > H - Delta;
    Tinterior = T > 0 & T <= H - Delta;   
    % Define the derivatives for the three regions
    ENGderivative = zeros(size(T));
    ENGderivative(Tzeros) = NaN;
    ENGderivative(Toob) = NaN;
    T = T(Tinterior);
    ENGderivative(Tinterior) = online.*(r.*(1 - 2.*pN).*mu0/2).*exp(-CTdiscrate.*T) - c.*r .*exp(-CTdiscrate.*T) ...
        + (1-pN).*(P_T(T).*sigmaZ(T).*normpdf(alphaN(T)./sigmaZ(T)) + P_T(T) .* (1 - normcdf(alphaN(T)./sigmaZ(T))).*mu0...
        + P(T).*sigmaZ_T(T).*normpdf(alphaN(T)./sigmaZ(T)) + P(T).*CTdiscrate.*(alphaN(T) + mu0)) ...
        + pN.*(P_T(T).*sigmaZ(T).*normpdf(alphaS(T)./sigmaZ(T)) - P_T(T).* (1 - normcdf(alphaS(T)./sigmaZ(T))).*mu0...
        + P(T).*sigmaZ_T(T).*normpdf(alphaS(T)./sigmaZ(T)) + P(T).*CTdiscrate.*(alphaS(T) - mu0));

end




