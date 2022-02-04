function [T,r,ENG] = TrasymptoticP(basic,P)
n0 = basic.n0;
sigmaX = basic.sigmaX;
mu0 = basic.mu0;
c = basic.c;
r = basic.r;
online = basic.online;
pN = basic.pN;
cr = basic.cr;

T = basic.Tmax;
constant = sqrt(n0.*sigmaX.^2) .* normpdf(sqrt(n0).*mu0./sigmaX)/...
    ((cr + c.*T - online.*(1 - 2.*pN).*mu0.*T./2).*T);
r = sqrt(constant.*P);
T = T.*ones(size(P));
ENG = P.*sigmaX./sqrt(n0).*(Psi(sqrt(n0).*mu0./sigmaX) + (1 - pN).*sqrt(n0).*mu0./sigmaX);

end

function y = Psi(x)
index = x == Inf;
y = zeros(size(x));
y(index) = 0;
x = x(~index);
y(~index) = normpdf(x) - x .* (1 - normcdf(x));
end