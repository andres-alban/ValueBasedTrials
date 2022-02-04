function [rstar,ENGstar,rdiscrete,ENGdiscrete] = OneShotMaxExpectedNetGainFreq(basic,globalsearch,Q)
% This function maximizes the the expected net gain when the adoption 
% decision is a two-sided frequentist hyposthesis test with respect to the
% recruitment rate when the sample size is constrained to be Q. If Q is not
% provided, the smallest sample size that guarantees errors and relevant
% difference specified in basic is used.
% Output:
% rstar and ENGstar are the maximizer and maximum of the continuous
% approximation using the built-in function fminbnd.
% rdiscrete and ENGdiscrete are the maximizer and maximum of the discrete
% problem, i.e. only values of of r are allowed such that rT=integer. These
% values are found through a simple search of the finite number of possible
% values of T. These values are only calculated if the output has four
% arguments. In that case, the rstar and rdiscrete are compared and a
% significant discrepancy gives a warning
%
% Revision: AA 9-4-2018
if nargin < 2
    globalsearch = false;
end
if nargin < 3
    ENG = @(r) OneShotExpectedNetGainFreq(basic,r); % Define the function to be maximized
else
    ENG = @(r) OneShotExpectedNetGainFreq(basic,r,Q); % Define the function to be maximized
end
ENGneg = @(r) -ENG(r);  % Minimize the negative instead


if ~globalsearch
    opts = optimset('TolX',min(basic.r./1e6,1e-6));
    [rstar,ENGstar,exitflag,output] = fminbnd(ENGneg,0,basic.zeta,opts);
    if exitflag ~= 1  % Output a warning if the optimization did not end properly
        warning('OneShotMaxExpectedNetGainFreq did not terminate properly')
        exitflag
        output
    end
else
    rstar = 0;
    ENGstar = ENGneg(rstar);
    for i = 1:5
        opts = optimoptions(@fmincon,'Algorithm','sqp','Display','off'); % Define the problem and options
        X0 = rand(1)*(basic.Tmax);
        [rstartemp,ENGstartemp,exitflag,output] = fmincon(ENGneg,X0,[],[],[],[],0,basic.Tmax,[],opts);
        if exitflag ~= 1  % Output a warning if the optimization did not end properly
            warning('OneShotMaxExpectedNetGainFreq did not terminate properly')
            exitflag
            output
        end
        if ENGstartemp < ENGstar
            rstar = rstartemp;
            ENGstar = ENGstartemp;
        end
    end
end
ENGstar = -ENGstar;  % Obtain the value of ENGstar


% Make sure that the maximum is larger than the value at at the
% endpoints
maxr = basic.zeta;
if ENGstar <= ENG(0)
    ENGstar = ENG(0);
    rstar = 0;
elseif ENGstar <= ENG(maxr)
    ENGstar = ENG(maxr);
    rstar = maxr;
end

% Calculate the discrete optimization
if nargout > 2
    h = basic.r/10;
    r = 0:h:(basic.zeta); % Define the grid of allowed values of r
    ENGvec = ENG(r); % Evaluate the function at the allowed values of r
    [ENGdiscrete,indexdiscrete] = max(ENGvec);  % Find the maximum
    rdiscrete = r(indexdiscrete);
    if abs(rstar - rdiscrete) > h  % If the difference between rdiscrete and rstar is larger than 1/r then there is an issue in the optimization
        warning('OneShotMaxExpectedNetGain does not match discrete result')
    end
end


end




