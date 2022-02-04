function [T,ENG] = TasymptoticH(basic,H)

n0 = basic.n0;
sigmaX = basic.sigmaX;
mu0 = basic.mu0;
c = basic.c;
r = basic.r;
online = basic.online;
pN = basic.pN;
zeta = basic.zeta;
Delta = basic.Delta;

if basic.CTdiscrate > 0
  warning('The approximation of TasymptoticH assumes zero discounting rate')
end

a1 = sqrt(n0).*mu0./sigmaX;
constant = c.*r - online.*r.*(1 - 2.*pN).*mu0/2 + zeta.*sqrt(sigmaX.^2./n0).*(Psi(a1) + (1 - pN).*a1);
constant = zeta.*sqrt(n0.*sigmaX.^2).*normpdf(a1)/(constant.*r);
T = sqrt(constant.*(H - Delta));

ENG = zeta.*(H - Delta).*sigmaX./sqrt(n0).*(Psi(sqrt(n0).*mu0./sigmaX) + (1 - pN).*sqrt(n0).*mu0./sigmaX);

end

function y = Psi(x)
index = x == Inf;
y = zeros(size(x));
y(index) = 0;
x = x(~index);
y(~index) = normpdf(x) - x .* (1 - normcdf(x));
end