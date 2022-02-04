function [T,ENG] = OneShotTasymptoticapprox(basic,PH)
% Compute the optimal T and expected net gain in the asymptote as P 
% (fixed patient pool) or H (fixed horizon) go to infinity, when r is fixed
% to basic.r.
  if basic.CTdiscrate > 0
    [T,ENG] = Tasymptoticdisc(basic)
    T = ones(size(PH)).*T;
    ENG = ones(size(PH)).*ENG;
  elseif strcmp(basic.horizon,'MktExcl') || strcmp(basic.horizon,'FixedPool')
    [T,ENG] = TasymptoticP(basic,PH);
  elseif strcmp(basic.horizon,'Patent') || strcmp(basic.horizon,'FixedHorizon')
    [T,ENG] = TasymptoticH(basic,PH);
  else
    error('No asymptotic approximation for your problem');
  end
end
