function [basic,rval] = MixInputModifier(basic,basicarray)
% This function modifies the model specified by basic by changing the value
% of the fields of structure basic. basicarray is a cell array that comes
% in pairs. The odd entries include a character vector specifying the field
% that has to be changed. The entries following the name of the field are
% the new value assigned to the field, e.g. a valid input is
% basicarray ={'c',2000,'horizon','MktExcl'}
%
% Revision: AA 28/6/2021


basicarrlen = length(basicarray);
rval = 1;
for i=1:(basicarrlen/2)  % For each pair
    if isfield(basic,basicarray{2*i-1})  % if the specified field exists
        basic.(basicarray{2*i-1}) = basicarray{2*i};  % change its value to the new assigned value
    else
        warning('invalid basic field: %s',char(basicarray{2*i-1}));
        rval = 0;
    end
end

% The function basic.P will be adjusted automatically based on the value of
% basic.horizon
if strcmp(basic.horizon,'Patent') || strcmp(basic.horizon,'FixedHorizon')
    basic.P = @(T) basic.zeta.*(basic.H-basic.Delta.*(T>0)-T);
    basic.P_T = @(T) -basic.zeta;
elseif strcmp(basic.horizon,'MktExcl') || strcmp(basic.horizon,'FixedPool')
    basic.P = @(T) basic.Population;
    basic.P_T = @(T) 0;
else
    warning('basic.horizon was specified to an unknown specification. Define manually basic.P and basic.P_T appropriately.')
end

% The function basic.ccap is adjusted if cfix or cr are changed
for i = 1:(basicarrlen/2)
   if strcmp(basicarray{2*i-1},'cfix') || strcmp(basicarray{2*i-1},'cr')
       basic.ccap = @(r) basic.cfix + basic.cr*r;
%        warning('basic.ccap has been assumed to be linear. Define manually basic.ccap for other desired functional forms')
       break
   end
end
