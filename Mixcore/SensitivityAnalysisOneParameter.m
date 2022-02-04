function [Tstar,ENGstar,PCDstar,p] = SensitivityAnalysisOneParameter(basic,parameter,parval)
% This is a generic function to perform a quick sensitivity analysis of the
% model with respect to a single parameter.
% Input:
% % basic: The model, defined usually with MixInputConstructor.m
% % paramter: a character vector with the name of the parameter
% % parval: the values of the parameter to be evaluated in the analysis
% Output:
% % Tstar: a vector containing the maximizer of the ENG when the parameter
%           takes the values in parval
% % ENGstar: a vector containing the maxima of ENG when the parameter
%           takes the values in parval
% % PCDstar: a vector containing the probability of correct decision at the
%           at the optimal trial length when the parameter takes the values
%           in parval
% % p = parval (this output is for consistency with
%               SensitivityAnalysisTwoParameters)
% 
% Revision: AA 25-11-2017

p = parval;
if isfield(basic,parameter) % Check that the parameter is in the model
    %Initialize the vectors
    Tstar = zeros(size(p));
    ENGstar = Tstar;
    PCDstar = Tstar;
    for i = 1:length(p)  % Loop over all the parameter values and find the maximum and maximizer of ENG
        basic = MixInputModifier(basic,{parameter,p(i)});
        [Tstar(i),ENGstar(i)] = OneShotMaxExpectedNetGain(basic);
        PCDstar(i) = OneShotPCD(basic,Tstar(i));
    end
else
    warning('invalid basic field: %s',parameter);  % Give a warning if the parameter does not exist
end
    