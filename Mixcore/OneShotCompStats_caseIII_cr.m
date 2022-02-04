function [dENG,dT,dr,dQ,Vrr] = OneShotCompStats_caseIII_cr(basic,T,r,cr,theta,h,hr)
% This function computes the numerical comparative statics for the setup 
% cost factor cr in case III, where the setup cost is given by cr*r^theta 
% for the problem defined by basic with optimal decision variables T and r.
% h, hT, and hr are the steps used to compute numerical derivatives.

V = OneShotExpectedNetGain(basic,T,r);
Vrr = (OneShotExpectedNetGain(basic,T,r+hr) - 2*V + OneShotExpectedNetGain(basic,T,r-hr))/(hr^2);
basic_hp = basic;
basic_hm = basic;
basic_hp.ccap = @(r) (cr+h)*r.^(theta);
basic_hm.ccap = @(r) (cr-h)*r.^(theta);
Vrparam = (OneShotExpectedNetGain(basic_hp,T,r+hr) - OneShotExpectedNetGain(basic_hm,T,r+hr) - OneShotExpectedNetGain(basic_hp,T,r-hr) + OneShotExpectedNetGain(basic_hm,T,r-hr))/(4*hr*h);
dENG = (OneShotExpectedNetGain(basic_hp,T,r) - OneShotExpectedNetGain(basic_hm,T,r))/(2*h);
dT = 0;
dr = -Vrparam/Vrr;
dQ = (dT*r+dr*T)/2;

end