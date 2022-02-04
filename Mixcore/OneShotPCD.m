function PCD = OneShotPCD(basic,T,r)
% This function calculates the probability of correct decision for a
% model specified by basic, trial length T, and recruitment rate r.

% At the moment the function does not work for unknown variance
if nargin < 3
    integrand = @(w) OneShotCPCD(w,basic,T) .* normpdf(w,basic.mu0,basic.sigmaX/sqrt(basic.n0));
else
    integrand = @(w) OneShotCPCD(w,basic,T,r) .* normpdf(w,basic.mu0,basic.sigmaX/sqrt(basic.n0));
end
PCD = integral(integrand,-inf,inf);



end

function PCD = OneShotPCDold(basic,T)
% Take parameters from the structure basic
c = basic.c;
pN = basic.pN;
mu0 = basic.mu0;
IS = basic.IS;
IN = basic.IN;
n0 = basic.n0;
sigmaX = basic.sigmaX;
r = basic.r;
P = basic.P;
CTdiscrate = basic.CTdiscrate;
zeta = basic.zeta;

% IN and IS required to have an adoption decision with type I error 0.05
% when the optimal trial length is 1.6
% IN = (1-pN).*P(1.6)./(n0+r.*1.6).*(sigmaX.*sqrt(r.*1.6).*norminv(1 - 0.025)+n0.*mu0);
% IS = pN.*P(1.6)./(n0+r.*1.6).*(sigmaX.*sqrt(r.*1.6).*norminv(1 - 0.025) - n0.*mu0);


% Define useful variables that are also defined in the paper
sigma0 = sqrt(sigmaX.^2./n0);

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
    alphaS = @(T) -IS ./ (pN .* Pdisc(T));
end

PNcorrect = zeros(size(T));
PScorrect = PNcorrect;
PMcorrect = PNcorrect;
PNincorrect = PNcorrect;
PSincorrect = PNcorrect;
PMincorrect = PNcorrect;
PCD = PNcorrect;
PID = PNcorrect;
sizeN = PNcorrect;
sizeS = PNcorrect;


for i = 1:length(T)
    if T(i) == 0
        if mu0 > alphaN(T(i))
            PCD(i) = 1 - normcdf((alphaN(T(i))-mu0)./sigma0);
        elseif mu0 < alphaS(T(i))
            PCD(i) = normcdf((alphaS(T(i))-mu0)./sigma0);
        else
            PCD(i) = normcdf((alphaN(T(i))-mu0)./sigma0) - normcdf((alphaS(T(i))-mu0)./sigma0);
        end
    else
        
        % Probability of adopting N given that N is the correct decision
        integrand = @(w) (1 - normcdf(sqrt(r.*T(i))./sigmaX.*((n0 + r.*T(i))./(r.*T(i)).*alphaN(T(i)) - n0.*mu0./(r.*T(i)) - w))).*normpdf((w - mu0)./sigma0);
        PNcorrect(i) = 1/(sigma0.*(1 - normcdf((alphaN(T(i)) - mu0)./sigma0))) .* integral(integrand,alphaN(T(i)),inf);
        
        % Probability of adopting S given that S is the correct decision
        integrand = @(w)  normcdf(sqrt(r.*T(i))./sigmaX.*((n0 + r.*T(i))./(r.*T(i)).*alphaS(T(i)) - n0.*mu0./(r.*T(i)) - w)) .* normpdf((w - mu0)./sigma0);
        PScorrect(i) = 1/(sigma0.* normcdf((alphaS(T(i)) - mu0)./sigma0)) .* integral(integrand,-inf,alphaS(T(i)));
        
        % Probability of adopting M given that M is the correct decision
        integrand = @(w) (normcdf(sqrt(r.*T(i))./sigmaX.*((n0 + r.*T(i))./(r.*T(i)).*alphaN(T(i)) - n0.*mu0./(r.*T(i)) - w)) ...
            - normcdf(sqrt(r.*T(i))./sigmaX.*((n0 + r.*T(i))./(r.*T(i)).*alphaS(T(i)) - n0.*mu0./(r.*T(i)) - w))) ...
            .* normpdf((w - mu0)./sigma0);
        PMcorrect(i) = 1/(sigma0.* (normcdf((alphaN(T(i)) - mu0)./sigma0) - normcdf((alphaS(T(i)) - mu0)./sigma0))) ...
            .* integral(integrand,alphaS(T(i)),alphaN(T(i)));
        
        % Probability of not adopting N given that N is the correct decision
        PNincorrect(i) = 1 - PNcorrect(i);
        
        % Probability of not adopting S given that S is the correct decision
        PSincorrect(i) =  1 - PScorrect(i);
        
        % Probability of not adopting M given that M is the correct decision
        PMincorrect(i) = 1 - PMcorrect(i);
        
        PCD(i) = (1 - normcdf((alphaN(T(i)) - mu0)./sigma0)).*PNcorrect(i) + normcdf((alphaS(T(i)) - mu0)./sigma0).*PScorrect(i)...
            + (normcdf((alphaN(T(i)) - mu0)./sigma0) - normcdf((alphaS(T(i)) - mu0)./sigma0)) .* PMcorrect(i);
        
        PID(i) = 1 - PCD(i);
        
        sizeN(i) = 1 - normcdf(sqrt(r.*T(i))./sigmaX.*((n0+r.*T(i)).*alphaN(T(i))-n0.*mu0)./(r.*T(i)));
        sizeS(i) = normcdf(sqrt(r.*T(i))./sigmaX.*((n0+r.*T(i)).*alphaS(T(i))-n0.*mu0)./(r.*T(i)));
    end
end

end

