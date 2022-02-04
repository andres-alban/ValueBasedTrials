function ProfherTimeHorizon(dosave)
close all
basic = ProfherSetParameters();
CTdiscrate = basic.CTdiscrate;
basic = MixInputModifier(basic,{'CTdiscrate',0});

H = [12,15,18,20:20:400];
P = basic.zeta.*H;

[Tstarmkt,ENGstarmkt,~,~] = SensitivityAnalysisOneParameter(basic,'Population',P);
Qstarmkt = Tstarmkt .* basic.r ./ 2;
[Tasymptmkt,ENGasymptmkt] = OneShotTasymptoticapprox(basic,P);
Qstarmktpred = 0.5*basic.r.*Tasymptmkt;
basic = MixInputModifier(basic,{'horizon','Patent'});
[Tstarpat,ENGstarpat,~,~] = SensitivityAnalysisOneParameter(basic,'H',H);
Qstarpat = Tstarpat .* basic.r ./ 2;
[Tasymptpat,ENGasymptpat] = OneShotTasymptoticapprox(basic,H);
Qstarpatpred = 0.5*basic.r.*Tasymptpat;




% same procedure with discounting
basic = MixInputModifier(basic,{'CTdiscrate',CTdiscrate,'horizon','MktExcl'});
[Tstarmktdisc,ENGstarmktdisc,~,~] = SensitivityAnalysisOneParameter(basic,'Population',P);
Qstarmktdisc = Tstarmktdisc .* basic.r ./ 2;
basic = MixInputModifier(basic,{'horizon','Patent'});
[Tstarpatdisc,ENGstarpatdisc,~,~] = SensitivityAnalysisOneParameter(basic,'H',H);
Qstarpatdisc = Tstarpatdisc .* basic.r ./ 2;
[Tasymptdisc,ENGasymptdisc] = OneShotTasymptoticapprox(basic,P);
Qstardiscpred = 0.5*basic.r.*Tasymptdisc;


figure(1)
plot(H./12,Qstarpat,'b-',H./12,Qstarpatdisc,'k--',H./12,Qstarpatpred,'r-.',H./12,Qstardiscpred,'r:','LineWidth',2)
xlabel('Time horizon H (years)');
ylabel('Optimal sample size Q^*_H')
ax = axis;
ax(4) = ax(4) + 5;
axis(ax);
legend('Undiscounted','Discounted','Asymptotic approx. undiscounted','Asymptotic approx. discounted','Location','SouthEast')
UtilStdizeFigureAbs(1);

figure(2)
plot(P,Qstarmkt,'b-',P,Qstarmktdisc,'k--',P,Qstarmktpred,'r-.',P,Qstardiscpred,'r:','LineWidth',2)
xlabel('P');
ylabel('Optimal sample size Q^*_P')
legend('Undiscounted','Discounted','Asymptotic approx. undiscounted','Asymptotic approx. discounted','Location','Northwest')
xlim([P(1) P(end)])
UtilStdizeFigureAbs(2);

figure(3)
plot(H./12,ENGstarpat./1e6,'b-',H./12,ENGstarpatdisc./1e6,'k--',H./12,ENGasymptpat./1e6,'r-.',H./12,ENGasymptdisc./1e6,'r:','LineWidth',2)
xlabel('Time horizon H (years)');
ylabel('Optimal expected net gain V^* (millions of £)')
legend('Undiscounted','Discounted','Asymptotic approx. undiscounted','Asymptotic approx. discounted','Location','NorthWest')
ax = axis;
ax(3) = 0;
ax(4) = ax(4) + 10;
axis(ax);
UtilStdizeFigureAbs(3);

figure(4)
plot(P,ENGstarmkt./1e6,'b-',P,ENGstarmktdisc./1e6,'k--',P,ENGasymptmkt./1e6,'r-.',P,ENGasymptdisc./1e6,'r:','LineWidth',2)
xlabel('P');
ylabel('Optimal expected net gain V^* (millions of £)')
legend('Undiscounted','Discounted','Asymptotic approx. undiscounted','Asymptotic approx. discounted','Location','SouthEast')
ax = axis;
ax(3) = 0;
axis(ax);
UtilStdizeFigureAbs(4);

% %%%%%%%%% CPCS plots
% H = [72 144];
% P = basic.zeta.*H;
% w = linspace(-3000,3000,3000);
% 
% basic = MixInputModifier(basic,{'horizon','MktExcl','Population',P(1)});
% [Tstar,~] = OneShotMaxExpectedNetGain(basic);
% CPCSmkt1 = OneShotCPCD(w,basic,Tstar);
% Powermkt1 = PowerPlotTwoSided(w,basic,Tstar,'Frequentist','Frequentist');
% 
% basic = MixInputModifier(basic,{'horizon','MktExcl','Population',P(2)});
% [Tstar,~] = OneShotMaxExpectedNetGain(basic);
% CPCSmkt2 = OneShotCPCD(w,basic,Tstar);
% Powermkt2 = PowerPlotTwoSided(w,basic,Tstar,'Frequentist','Frequentist');
% 
% basic = MixInputModifier(basic,{'horizon','Patent','H',H(1)});
% [Tstar,~] = OneShotMaxExpectedNetGain(basic);
% CPCSpat1 = OneShotCPCD(w,basic,Tstar);
% Powerpat1 = PowerPlotTwoSided(w,basic,Tstar,'Frequentist','Frequentist');
% 
% basic = MixInputModifier(basic,{'horizon','Patent','H',H(2)});
% [Tstar,~] = OneShotMaxExpectedNetGain(basic);
% CPCSpat2 = OneShotCPCD(w,basic,Tstar);
% Powerpat2 = PowerPlotTwoSided(w,basic,Tstar,'Frequentist','Frequentist');
% 
% figure(4)
% plot(w,CPCSpat2 - CPCSpat1,'b-',w,CPCSmkt2 - CPCSmkt1,'r--','LineWidth',2)
% delta = basic.relevantdifference;
% hold on
% plot([delta delta],[0 0.06],'k--',[-delta -delta],[0 0.06],'k--','LineWidth',1)
% hold off
% legend('\Delta CPCS_{H}','\Delta CPCS_{P}','Location','NorthEast')
% xlabel('w')
% ylabel('CPCS')
% UtilStdizeFigureAbs(4);

% Save the figures if dosave=true
if nargin < 1
    dosave = false;
end
if dosave
    UtilSaveFigFile(1, 'Profher', 'ProfherTimeHorizon_pat',[],'eps');
    UtilSaveFigFile(2, 'Profher', 'ProfherTimeHorizon_mkt',[],'eps');
    UtilSaveFigFile(3, 'Profher', 'ProfherTimeHorizon_ENG_pat',[],'eps');
    UtilSaveFigFile(4, 'Profher', 'ProfherTimehorizon_ENG_mkt',[],'eps');
end








