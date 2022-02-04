function [Designstar,ENGstar,Designdiscrete,ENGdiscrete] = OneShotMaxExpectedNetGain(basic,globalsearch,optim_variables,iters,verbose)
% This function maximizes the the expected net gain with respect to the
% trial length, recruitment rate or both.
% Output:
% Designstar and ENGstar are the maximizer and maximum of the continuous
% approximation using the built-in function fmincon.
% Designdiscrete and ENGdiscrete are the maximizer and maximum of the discrete
% problem, i.e. only values of of T are allowed such that rT=integer. These
% values are found through a simple search of the finite number of possible
% values of T. These values are only calculated if the output has four
% arguments. In that case, the Tstar and Tdiscrete are compared and a
% significant discrepancy gives a warning
%
% Revision: AA 28-6-2021
if nargin < 2
    globalsearch = false;
    optim_variables = 'T';
elseif nargin < 3
    optim_variables = 'T';
end
if optim_variables == true % for backwards compaibility
    optim_variables = 'Tr';
elseif optim_variables == false
    optim_variables = 'T';
end
if nargin < 4
    iters = 100;
end
if nargin < 5
    verbose = true;
end

if strcmp(optim_variables,'T') % optimize over T only
    ENG = @(T) OneShotExpectedNetGain(basic,T); % Define the function to be maximized
    ENGneg = @(T) -ENG(T);  % Minimize the negative instead
    
    
    if ~globalsearch
        opts = optimset('TolX',min(1/basic.r,1e-6),'Display','off');
        [Designstar,ENGstar,exitflag,output] = fminbnd(ENGneg,0,basic.Tmax,opts);
        if exitflag ~= 1 && verbose  % Output a warning if the optimization did not end properly
            warning('OneShotMaxExpectedNetGain did not terminate properly')
            exitflag
            output
        end
    else
        Designstar = 0;
        ENGstar = ENGneg(Designstar);
        for i = 1:iters
            opts = optimoptions(@fmincon,'Algorithm','sqp','Display','off'); % Define the problem and options
            X0 = rand(1)*(basic.Tmax);
            [Tstartemp,ENGstartemp,exitflag,output] = fmincon(ENGneg,X0,[],[],[],[],0,basic.Tmax,[],opts);
            if exitflag ~= 1 && verbose  % Output a warning if the optimization did not end properly
                warning('OneShotMaxExpectedNetGain did not terminate properly')
                exitflag
                output
            end
            if ENGstartemp < ENGstar
                Designstar = Tstartemp;
                ENGstar = ENGstartemp;
            end
        end
    end
    ENGstar = -ENGstar;  % Obtain the value of ENGstar
    
    
    % Make sure that the maximum is larger than the value at at the
    % endpoints
    maxT = basic.Tmax;
    if ENGstar <= ENG(0)
        ENGstar = ENG(0);
        Designstar = 0;
    end
    if ENGstar <= ENG(maxT)
        ENGstar = ENG(maxT);
        Designstar = maxT;
    end
    
    % Calculate the discrete optimization
    if nargout > 2
        h = 1/basic.r;
        T = 0:h:(basic.Tmax); % Define the grid of allowed values of T
        ENGvec = ENG(T); % Evaluate the function at the allowed values of T
        [ENGdiscrete,indexdiscrete] = max(ENGvec);  % Find the maximum
        Designdiscrete = T(indexdiscrete);
        if abs(Designstar - Designdiscrete) > h  % If the difference between Tdiscrete and Tstar is larger than 1/r then there is an issue in the optimization
            warning('OneShotMaxExpectedNetGain does not match discrete result')
        end
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(optim_variables,'r') % optimize over r only
    ENG = @(r) OneShotExpectedNetGain(basic,basic.T,r); % Define the function to be maximized
    ENGneg = @(r) -ENG(r);  % Minimize the negative instead
    
    
    if ~globalsearch
        opts = optimset('TolX',min(1/basic.T,1e-6),'Display','off');
        [Designstar,ENGstar,exitflag,output] = fminbnd(ENGneg,0,basic.rmax,opts);
        if exitflag ~= 1 && verbose  % Output a warning if the optimization did not end properly
            warning('OneShotMaxExpectedNetGain did not terminate properly')
            exitflag
            output
        end
    else
        Designstar = 0;
        ENGstar = ENGneg(Designstar);
        for i = 1:iters
            opts = optimoptions(@fmincon,'Algorithm','sqp','Display','off'); % Define the problem and options
            X0 = rand(1)*(basic.rmax);
            [rstartemp,ENGstartemp,exitflag,output] = fmincon(ENGneg,X0,[],[],[],[],0,basic.rmax,[],opts);
            if exitflag ~= 1 && verbose  % Output a warning if the optimization did not end properly
                warning('OneShotMaxExpectedNetGain did not terminate properly')
                exitflag
                output
            end
            if ENGstartemp < ENGstar
                Designstar = rstartemp;
                ENGstar = ENGstartemp;
            end
        end
    end
    ENGstar = -ENGstar;  % Obtain the value of ENGstar
    
    
    % Make sure that the maximum is larger than the value at at the
    % endpoints
    maxr = basic.rmax;
    if ENGstar <= ENG(0)
        ENGstar = ENG(0);
        Designstar = 0;
    end
    if ENGstar <= ENG(maxr)
        ENGstar = ENG(maxr);
        Designstar = maxr;
    end
    
    % Calculate the discrete optimization
    if nargout > 2
        h = 1/basic.T;
        r = 0:h:(basic.rmax); % Define the grid of allowed values of T
        ENGvec = zeros(size(r));
        for i =1:length(r)
            ENGvec(i) = ENG(r(i)); % Evaluate the function at the allowed values of T
        end
        [ENGdiscrete,indexdiscrete] = max(ENGvec);  % Find the maximum
        Designdiscrete = r(indexdiscrete);
        if abs(Designstar - Designdiscrete) > h  % If the difference between Tdiscrete and Tstar is larger than 1/r then there is an issue in the optimization
            warning('OneShotMaxExpectedNetGain does not match discrete result')
        end
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

else  % optimize over both T and r
    ENG = @(X) OneShotExpectedNetGain(basic,X(1),X(2)); % Define the function to be maximized
    ENGneg = @(X) -ENG(X);  % Minimize the negative instead
    
    
    if ~globalsearch
        opts = optimoptions(@fmincon,'Algorithm','sqp','Display','off'); % Define the problem and options
        X0 = rand(2,1).*[basic.Tmax;basic.rmax];
        [Designstar,ENGstar,exitflag,output] = fmincon(ENGneg,X0,[],[],[],[],[0;0],[basic.Tmax;basic.rmax],[],opts);
        if exitflag < 1 && verbose  % Output a warning if the optimization did not end properly
            warning('OneShotMaxExpectedNetGain did not terminate properly')
            exitflag
            output
        end
    else
        Designstar = [0;0];
        ENGstar = ENGneg(Designstar);
        for i = 1:iters
            opts = optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',10,'Display','off'); % Define the problem options
            X0 = rand(2,1).*[basic.Tmax;basic.rmax];
            [Tstartemp,ENGstartemp] = fmincon(ENGneg,X0,[],[],[],[],[0;0],[basic.Tmax;basic.rmax],[],opts);
            if ENGstartemp < ENGstar
                Designstar = Tstartemp;
                ENGstar = ENGstartemp;
            end
        end
        opts = optimoptions(@fmincon,'Algorithm','sqp','Display','off'); % Define the problem options
        [Designstar,ENGstar,exitflag,output] = fmincon(ENGneg,Designstar,[],[],[],[],[0;0],[basic.Tmax;basic.rmax],[],opts);
        if exitflag < 1 && verbose  % Output a warning if the optimization did not end properly
            warning('OneShotMaxExpectedNetGain did not terminate properly')
            exitflag
            output
        end
    end
    ENGstar = -ENGstar;  % Obtain the value of ENGstar
    
    
    % Make sure that the maximum is larger than the value at at the
    % endpoints
    maxT = basic.Tmax;
    maxr = basic.rmax;
    if ENGstar <= ENG([0;0])
        ENGstar = ENG([0;0]);
        Designstar = [0;0];
    end
    if ENGstar <= ENG([maxT;maxr])
        ENGstar = ENG([maxT;maxr]);
        Designstar = [maxT;maxr];
    end
    
    % Calculate the discrete optimization
    if nargout > 2
        warning('The discrete optimization does not work when optimizing over T and r.')
        ENGdiscrete = [];
        Designdiscrete = [];
    end
    
end




