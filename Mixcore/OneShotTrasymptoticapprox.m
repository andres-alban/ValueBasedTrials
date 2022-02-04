function [T,r,ENG] = OneShotTrasymptoticapprox(basic,PH)
% Compute the optimal T, r, and expected net gain in the asymptote as P 
% (fixed patient pool) or H (fixed horizon) go to infinity.
  if basic.CTdiscrate > 0
    [T,r,ENG] = Trasymptoticdisc(basic);
    T = ones(size(PH)).*T;
    r = ones(size(PH)).*r;
    ENG = ones(size(PH)).*ENG;
  elseif strcmp(basic.horizon,'MktExcl') || strcmp(basic.horizon,'FixedPool')
    [T,r,ENG] = TrasymptoticP(basic,PH);
  elseif strcmp(basic.horizon,'Patent') || strcmp(basic.horizon,'FixedHorizon')
    [T,r,ENG] = TrasymptoticH(basic,PH);
  else
    error('No asymptotic approximation for your problem');
  end
end
