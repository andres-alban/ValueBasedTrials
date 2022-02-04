function [Tstar,ENGstar,PCDstar,p1,p2] = SensitivityAnalysisTwoParameters(basic,parameter1,parval1,parameter2,parval2)
% This is a generic function to perform a quick sensitivity analysis of the
% model with respect to two parameters.
% Input:
% % basic: The model, defined usually with MixInputConstructor.m
% % paramter1: a character vector with the name of the first parameter
% % parval1: the values of parameter1 to be evaluated in the analysis
% % paramter2: a character vector with the name of the second parameter
% % parval2: the values of parameter2 to be evaluated in the analysis
% Output:
% % Tstar: a matrix containing the maximizer of the ENG when the parameter
%           take the values in parval1 and parval2
% % ENGstar: a matrix containing the maxima of ENG when the parameters
%           take the values in parval1 and parval2
% % PCDstar: a matrix containing the probability of correct decision at the
%           at the optimal trial length when the parameter takes the values
%           in parval
% % p1: matrix of the values taken by parval1 for 3D plotting
% % p2: matrix of the values taken by parval2 for 3D plotting
% 
% The entries of Tstar and ENGstar match the solutions to the model with
% the parameters specified by the entries in p1 and p2. This can be easily
% plotted in 3D with commands such as surf(p1,p2,Tstar)
% 
% Revision: AA 25-11-2017


if isfield(basic,parameter1) && isfield(basic,parameter2) % check that the parameters exist
    [p1,p2] = meshgrid(parval1,parval2);  % create the matrices of parameter values
%     Initialize matrices of solutions
    Tstar = zeros(size(p1));
    ENGstar = Tstar;
    PCDstar = Tstar;
    [row,col] = size(p1);
    % iterate over all combinations of parameter values and find the
    % maximizers and maxima
    for i = 1:row
        for j = 1:col
        basic = MixInputModifier(basic,{parameter1,p1(i,j),parameter2,p2(i,j)});
        [Tstar(i,j),ENGstar(i,j)] = OneShotMaxExpectedNetGain(basic);
        PCDstar(i,j) = OneShotPCD(basic,Tstar(i,j));
        end
    end
else
    warning('invalid basic field: %s or %s',parameter1,parameter2);  % give a warning in case of an inexisting parameter
end
    