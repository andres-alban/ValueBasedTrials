function basic = MixInputConstructor(basicarray)
% This function generates the structure basic that defines the parameters
% of the model. It is a required input for most of other functions. The
% default output creates the model used in the Profher Paper. The model can
% be modified by inputting basicarray, a cell array that includes the name
% of the field, and the new value assigned to that field. This is done
% through the MixInputModifier.m function.
%
% Revision: AA 18-12-2017



% Parameters
basic.c = 1000;  % Marginal cost of an additional patient in the trial
basic.pN = 0.39;  % Proportion of patient in the new treatment with mix treatment
basic.mu0 = 0;  %  Mean of prior distribution for uknown mean reward
basic.IS = 0;  % Switching cost to standard treatment
basic.IN = 10^7;  % Switching cost to new treatment
basic.zeta = 7000;  % incidence rate
basic.Tmax = 20;
basic.rmax = 7000;
basic.H = 20;  % Time horizon
basic.Delta = 1;  % Delay to acquire info
basic.Population = 140000;  % Population to be treated after an adoption decision
basic.n0 = 20;  % Effective number of samples
basic.sigmaX = 2100;  % Known std deviation of each observation
basic.cfix = 110000;  % Fixed cost of trial
basic.cr = 42000; % Marginal cost of an additional unit of recruitment rate
basic.cT = 0; % Marginal cost of an additional unit of time in the trial
basic.r = 50;  % Recruitment rate
basic.T = 2; % trial duration
basic.horizon = 'FixedPool'; % 'FixedHorizon' || 'Patent' -> Model fixed horizon || Patent protection; 'FixedPool' || 'MktExcl' -> Model of fixed patient pool || market exclusivity
basic.online = 1; % 1 -> online | 0 -> offline
basic.CTdiscrate = 0; % Discount rate. If 0, undiscounted 
basic.relevantdifference = 500;  % minimum clinically relevant difference
basic.alpha = 0.05;  % Type I error specified for frequentist sample size calculations
basic.power = 0.8;  % Pre-specified power of the trial
basic.UnkVariance = false; % false if variance is known
basic.UnkVarianceShape = 3; % Shape parameter of inverse gamma distribution of unknown variance 
% basic.chi0 = 441000; % Scale parameter of inverse gamma distribution of unknown variance

%These functions are automatically updated based on the value of basic.horizon through MixInputModifier.m
basic.P = @(T) basic.zeta.*(basic.H - (T>0) * basic.Delta - T);  % Number of patients to be treated as a function of trial length. 
basic.P_T = @(T) -basic.zeta; % Derivative of P(T) in the interior
basic.ccap = @(r) basic.cfix + basic.cr*r;



if nargin == 1 && ~isempty(basicarray)
    basic = MixInputModifier(basic,basicarray); % Modify the model based on the input
end
