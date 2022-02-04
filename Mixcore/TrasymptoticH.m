function [T,r,ENG] = TrasymptoticH(basic,H)
n0 = basic.n0;
sigmaX = basic.sigmaX;
mu0 = basic.mu0;
c = basic.c;
online = basic.online;
pN = basic.pN;
zeta = basic.zeta;
Delta = basic.Delta;
cr = basic.cr;

if (c - online*(1 - 2*pN)*mu0/2) > 0
    K = sigmaX./sqrt(n0).*zeta.*(Psi(sqrt(n0).*mu0./sigmaX) + (1 - pN).*sqrt(n0).*mu0./sigmaX)./cr;
    constant = sqrt(n0).*sigmaX .* normpdf(sqrt(n0).*mu0./sigmaX)./...
        (c - online.*(1 - 2.*pN).*mu0/2);
    T = (constant./K.^2.*zeta.*(H - Delta)).^(1./4);
    r = (constant.*K.^2.*zeta.*(H - Delta)).^(1./4);
elseif (c - online*(1 - 2*pN)*mu0/2) == 0
    K = sigmaX./sqrt(n0).*zeta.*(Psi(sqrt(n0).*mu0./sigmaX) + (1 - pN).*sqrt(n0).*mu0./sigmaX)./cr;
    constant = sqrt(n0).*sigmaX .* normpdf(sqrt(n0).*mu0./sigmaX)./cr;
    T = (constant./K.^2.*zeta.*(H - Delta)).^(1/3);
    r = T.*K;
else
    T = H-Delta;
    r = zeta.*onses(size(T));
end

ENG = zeta.*(H - Delta).*sigmaX./sqrt(n0).*(Psi(sqrt(n0).*mu0./sigmaX) + (1 - pN).*sqrt(n0).*mu0./sigmaX);

end

function y = Psi(x)
index = x == Inf;
y = zeros(size(x));
y(index) = 0;
x = x(~index);
y(~index) = normpdf(x) - x .* (1 - normcdf(x));
end