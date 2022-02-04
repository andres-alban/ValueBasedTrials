function [T,ENG] = TasymptoticP(basic,P)
n0 = basic.n0;
sigmaX = basic.sigmaX;
mu0 = basic.mu0;
c = basic.c;
r = basic.r;
online = basic.online;
pN = basic.pN;

if basic.CTdiscrate > 0
  warning('The approximation of TasymptoticP assumes zero discounting rate')
end

constant = sqrt(n0.*sigmaX.^2) .* normpdf(sqrt(n0).*mu0./sigmaX)/...
    (r.^2.*(c - online.*(1 - 2.*pN).*mu0./2));
T = sqrt(constant.*P);
ENG = P.*sigmaX./sqrt(n0).*(Psi(sqrt(n0).*mu0./sigmaX) + (1 - pN).*sqrt(n0).*mu0./sigmaX);

end

function y = Psi(x)
index = x == Inf;
y = zeros(size(x));
y(index) = 0;
x = x(~index);
y(~index) = normpdf(x) - x .* (1 - normcdf(x));
end