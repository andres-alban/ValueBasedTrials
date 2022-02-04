function [dENG,dT,dr,dQ,Vrr] = OneShotCompStats_caseIII(basic,T,r,param,h,hr)
% This function computes the numerical comparative statics for parameter 
% param in case III for the problem defined by basic with optimal 
% decision variables T and r. h and hr are the steps used to compute
% numerical derivatives.

V = OneShotExpectedNetGain(basic,T,r);
Vrr = (OneShotExpectedNetGain(basic,T,r+hr) - 2*V + OneShotExpectedNetGain(basic,T,r-hr))/(hr^2);
basic_hp = MixInputModifier(basic,{param,basic.(param)+h});
basic_hm = MixInputModifier(basic,{param,basic.(param)-h});
Vrparam = (OneShotExpectedNetGain(basic_hp,T,r+hr) - OneShotExpectedNetGain(basic_hm,T,r+hr) - OneShotExpectedNetGain(basic_hp,T,r-hr) + OneShotExpectedNetGain(basic_hm,T,r-hr))/(4*hr*h);
dENG = (OneShotExpectedNetGain(basic_hp,T,r) - OneShotExpectedNetGain(basic_hm,T,r))/(2*h);
dT = 0;
dr = -Vrparam/Vrr;
dQ = (dT*r+dr*T)/2;

end