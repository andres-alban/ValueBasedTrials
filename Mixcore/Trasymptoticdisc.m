function [T,r,ENG] = Trasymptoticdisc(basic)
if basic.CTdiscrate == 0
    error('The input to Tasymptoticdisc has to be a discounted problem (CTdiscrate>0)');
end
basic = MixInputModifier(basic,{'horizon','MktExcl','Population',Inf});
[Tr,ENG] = OneShotMaxExpectedNetGain(basic,true,true);
T = Tr(1);
r = Tr(2);

