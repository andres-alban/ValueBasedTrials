function [T,ENG] = Tasymptoticdisc(basic)
if basic.CTdiscrate == 0
    error('The input to Tasymptoticdisc has to be a discounted problem (CTdiscrate>0)');
end
basic = MixInputModifier(basic,{'horizon','MktExcl','Population',Inf});
[T,ENG] = OneShotMaxExpectedNetGain(basic);



% n0 = basic.n0;
% sigmaX = basic.sigmaX;
% mu0 = basic.mu0;
% c = basic.c;
% r = basic.r;
% online = basic.online;
% pN = basic.pN;
% zeta = basic.zeta;
% IN = basic.IN;
% IS = basic.IS;
% rho = basic.CTdiscrate;
% Delta = basic.Delta;

% alphaN = IN./(zeta./rho.*(1 - pN));
% alphaS = IS./(zeta./rho.*pN);
% sigmaZ = @(T) sqrt(sigmaX^2 .* r .* T ./ (n0 .* (2.*n0 + r .* T)));
% sigmaZ_T = @(T) sqrt(n0.*sigmaX^2.*r ./ ((2.*n0 + r.*T).^3 .* T)); % Derivative of sigmaZ with respect to T
% 
% K = @(T) exp(-rho.*Delta)./(r.*(c - online.*(1 - 2.*pN).*mu0./2))...
%     .* ((1 - pN).* (sigmaZ_T(T).*normpdf((alphaN - mu0)./sigmaZ(T)) - sigmaZ(T).*rho.*Psi((alphaN-mu0)./sigmaZ(T)))...
%     + pN.* (sigmaZ_T(T).*normpdf((alphaS + mu0)./sigmaZ(T)) - sigmaZ(T).*rho.*Psi((alphaS + mu0)./sigmaZ(T))) );
% K2 = @(T) 1 - rho./(zeta.*K(T));
% K3 = @(T) K(T) - rho./zeta;
% TL = 1.98.*n0./(0.01.*r);
% TH = 0.02.*n0./(0.99.*r);
% Tasympt = fzero(K3,[TL,TH]);

% Tzeros = T == 0;
% Toob = T < 0;
% Tinterior = T > 0;% & T <= Tasympt;
% Tlarge = T > Tasympt;
% P = zeros(size(T));
% P(Tzeros) = 0;
% P(Toob) = NaN;
% P(Tinterior) = -zeta./rho.*log(K2(T(Tinterior)));
% P(Tlarge) = Tasympt;

% end
% 
% 
% function y = Psi(x)
% index = x == Inf;
% y = zeros(size(x));
% y(index) = 0;
% x = x(~index);
% y(~index) = normpdf(x) - x .* (1 - normcdf(x));
% end