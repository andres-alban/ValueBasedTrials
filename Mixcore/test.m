% This script just tests the use of some of the functions. It can be also
% useful to understand the way some of the functions are implemented.
basic = MixInputConstructor({'CTdiscrate',0.035});
ENG = @(T) OneShotExpectedNetGain(basic,T);
ENGderivative = @(T) OneShotExpectedNetGainDeriv(basic,T);
[Tstar,ENGstar,Tdiscrete,ENGdiscrete] = OneShotMaxExpectedNetGain(basic);
figure(1)
T = 0:(1/basic.r):(basic.H-basic.Delta);
plot(T,ENG(T),'.',Tstar,ENGstar,'ko',Tdiscrete,ENGdiscrete,'rx')
figure(2)
plot(T,ENGderivative(T),Tstar,ENGderivative(Tstar),'ko')

pN = 0.1:0.1:0.9;
[Tstar1,ENGstar1,pN] = SensitivityAnalysisOneParameter(basic,'pN',pN);
figure(3)
subplot(1,2,1)
plot(pN,Tstar1)
xlabel('p_N');ylabel('T^*')
subplot(1,2,2)
plot(pN,ENGstar1)
xlabel('p_N');ylabel('Maximum expected net gain')

[Tstar2,ENGstar2,p1,p2] = SensitivityAnalysisTwoParameters(basic,'mu0',-100:10:100,'n0',1:20:200);
figure(4)
subplot(1,2,1)
surf(p1,p2,Tstar2)
xlabel('mu0');ylabel('n0');zlabel('T^*')
subplot(1,2,2)
surf(p1,p2,ENGstar2)
xlabel('mu0');ylabel('n0');zlabel('Maximum expected net gain')