function Profherpn(dosave)
close all

basic = ProfherSetParameters();
I = [10e6 100e6];
pN = 0.01:0.01:0.5;


basic = MixInputModifier(basic,{'IN',I(1),'IS',I(1)});
[Tstar,ENGstarmkt1,~,~] = SensitivityAnalysisOneParameter(basic,'pN',pN);
Qstarmkt1 = basic.r .* Tstar;

basic = MixInputModifier(basic,{'horizon','Patent'});
[Tstar,ENGstarpat1,~,~] = SensitivityAnalysisOneParameter(basic,'pN',pN);
Qstarpat1 = basic.r .* Tstar;

basic = MixInputModifier(basic,{'horizon','MktExcl','IN',I(2),'IS',I(2)});
[Tstar,ENGstarmkt2,~,~] = SensitivityAnalysisOneParameter(basic,'pN',pN);
Qstarmkt2 = basic.r .* Tstar;

basic = MixInputModifier(basic,{'horizon','Patent','IN',I(2),'IS',I(2)});
[Tstar,ENGstarpat2,~,~] = SensitivityAnalysisOneParameter(basic,'pN',pN);
Qstarpat2 = basic.r .* Tstar;

figure(1)
plot(pN,Qstarpat1,'b-',pN,Qstarpat2,'r--',pN,Qstarmkt1,'k:',pN,Qstarmkt2,'g-.','LineWidth',2)
xlabel('Fraction in new treatment (p_N)');ylabel('Optimal number of pairwise allocations Q^*')
lgd = legend(['Q^*_{H}, I_N=I_S=',num2str(I(1)/1e6),'Mill.'],['Q^*_{H}, I_N=I_S=',num2str(I(2)/1e6),'Mill.'],...
    ['Q^*_{P}, I_N=I_S=',num2str(I(1)/1e6),'Mill.'],['Q^*_{P}, I_N=I_S=',num2str(I(2)/1e6),'Mill.']);
set(lgd,'Position',[0.558333339974994 0.230295570123763 0.295714279072625 0.252052537913393]);
UtilStdizeFigureAbs(1);

figure(2)
plot(pN,ENGstarpat1,'b-',pN,ENGstarpat2,'r--',pN,ENGstarmkt1,'k:',pN,ENGstarmkt2,'g-.','LineWidth',2)
xlabel('Fraction in new treatment (p_N)');ylabel('Optimal expected net gain V^*')
lgd = legend(['V^*_{H}, I_N=I_S=',num2str(I(1)/1e6),'Mill.'],['V^*_{H}, I_N=I_S=',num2str(I(2)/1e6),'Mill.'],...
    ['V^*_{P}, I_N=I_S=',num2str(I(1)/1e6),'Mill.'],['V^*_{P}, I_N=I_S=',num2str(I(2)/1e6),'Mill.']);
set(lgd,'Position',[0.570660646387833 0.276583514099783 0.291349809885932 0.217190889370933]);
UtilStdizeFigureAbs(2);

% Save the figures if dosave=true
if nargin < 1
    dosave = false;
end
if dosave
    UtilSaveFigFile(1, 'Profher', 'Profherpn_Q',[],'eps');
    UtilSaveFigFile(2, 'Profher', 'Profherpn_B',[],'eps');
end



