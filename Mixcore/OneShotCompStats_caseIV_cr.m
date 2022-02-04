function [dENG,dT,dr,dQ,H] = OneShotCompStats_caseIV_cr(basic,T,r,cr,theta,h,hT,hr)
% This function computes the numerical comparative statics for the setup 
% cost factor cr in case IV, where the setup cost is given by cr*r^theta 
% for the problem defined by basic with optimal decision variables T and r.
% h, hT, and hr are the steps used to compute numerical derivatives.

V = OneShotExpectedNetGain(basic,T,r);
VTT = (OneShotExpectedNetGain(basic,T+hT,r) - 2*V + OneShotExpectedNetGain(basic,T-hT,r))/(hT^2);
Vrr = (OneShotExpectedNetGain(basic,T,r+hr) - 2*V - OneShotExpectedNetGain(basic,T,r-hr))/(hr^2);
VTr = (OneShotExpectedNetGain(basic,T+hT,r+hr) - OneShotExpectedNetGain(basic,T+hT,r-hr) - OneShotExpectedNetGain(basic,T-hT,r+hr) + OneShotExpectedNetGain(basic,T-hT,r-hr))/(4*hT*hr);
H = det([VTT,VTr;VTr,Vrr]);
basic_hp = basic;
basic_hm = basic;
basic_hp.ccap = @(r) (cr+h)*r.^(theta);
basic_hm.ccap = @(r) (cr-h)*r.^(theta);
VTparam = (OneShotExpectedNetGain(basic_hp,T+hT,r) - OneShotExpectedNetGain(basic_hm,T+hT,r) - OneShotExpectedNetGain(basic_hp,T-hT,r) + OneShotExpectedNetGain(basic_hm,T-hT,r))/(4*hT*h);
Vrparam = (OneShotExpectedNetGain(basic_hp,T,r+hr) - OneShotExpectedNetGain(basic_hm,T,r+hr) - OneShotExpectedNetGain(basic_hp,T,r-hr) + OneShotExpectedNetGain(basic_hm,T,r-hr))/(4*hr*h);
dENG = (OneShotExpectedNetGain(basic_hp,T,r) - OneShotExpectedNetGain(basic_hm,T,r))/(2*h);
dT = -det([VTparam,VTr;Vrparam,Vrr])/H;
dr = -det([VTT,VTparam;VTr,Vrparam])/H;
dQ = (dT*r+dr*T)/2;

end