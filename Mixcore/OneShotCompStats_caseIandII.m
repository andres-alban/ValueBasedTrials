function [dENG,dT,dr,dQ,VTT] = OneShotCompStats_caseIandII(basic,T,r,param,h,hT)
% This function computes the numerical comparative statics for parameter 
% param in cases I and II for the problem defined by basic with optimal 
% decision variables T and r. h and hT are the steps used to compute
% numerical derivatives.

V = OneShotExpectedNetGain(basic,T,r);
VTT = (OneShotExpectedNetGain(basic,T+hT,r) - 2*V + OneShotExpectedNetGain(basic,T-hT,r))/(hT^2);
basic_hp = MixInputModifier(basic,{param,basic.(param)+h});
basic_hm = MixInputModifier(basic,{param,basic.(param)-h});
VTparam = (OneShotExpectedNetGain(basic_hp,T+hT,r) - OneShotExpectedNetGain(basic_hm,T+hT,r) - OneShotExpectedNetGain(basic_hp,T-hT,r) + OneShotExpectedNetGain(basic_hm,T-hT,r))/(4*hT*h);
dENG = (OneShotExpectedNetGain(basic_hp,T,r) - OneShotExpectedNetGain(basic_hm,T,r))/(2*h);
dT = -VTparam/VTT;
dr = 0;
dQ = (dT*r+dr*T)/2;

end