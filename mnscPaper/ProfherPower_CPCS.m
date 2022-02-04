function ProfherPower_CPCS(dosave)
close all
basic = ProfherSetParameters();
[Tstar,~] = OneShotMaxExpectedNetGain(basic);
delta = basic.relevantdifference;
w = linspace(-1.5*delta,1.5*delta,delta/2);
CPCSmkt = OneShotCPCD(w,basic,Tstar);
Powermkt = PowerPlotTwoSided(w,basic,Tstar,[],'Frequentist','Frequentist');

basic = MixInputModifier(basic,{'horizon','Patent'});
[Tstar,~] = OneShotMaxExpectedNetGain(basic);
CPCSpat = OneShotCPCD(w,basic,Tstar);
Powerpat = PowerPlotTwoSided(w,basic,Tstar,[],'Frequentist','Frequentist');

basic = MixInputModifier(basic,{'horizon','MktExcl','Population',2*basic.Population});
[Tstar2,~] = OneShotMaxExpectedNetGain(basic);
CPCSmkt2 = OneShotCPCD(w,basic,Tstar2);
Powermkt2 = PowerPlotTwoSided(w,basic,Tstar2,[],'Frequentist','Frequentist');

basic = MixInputModifier(basic,{'horizon','Patent','H',2*basic.H});
[Tstar2,~] = OneShotMaxExpectedNetGain(basic);
CPCSpat2 = OneShotCPCD(w,basic,Tstar2);
Powerpat2 = PowerPlotTwoSided(w,basic,Tstar2,[],'Frequentist','Frequentist');

figure(1)
plot(w,CPCSpat,'b-',w,CPCSpat2,'b-.',w,CPCSmkt,'r--',w,CPCSmkt2,'r:','LineWidth',2)
xlim([-1.5*delta,1.5*delta])
ylim([0,1])
hold on
plot([delta delta],[0 1],'k--',[-delta -delta],[0 1],'k--','LineWidth',1)
hold off
legend('CPCS_{H}','CPCS_{H} 2\times H','CPCS_{P}','CPCS_{P} 2\times P','Location','South')
xlabel('w')
ylabel('CPCS')
UtilStdizeFigureAbs(1);

figure(2)
plot(w,Powerpat,'b-',w,Powerpat2,'b-.',w,Powermkt,'r--',w,Powermkt2,'r:','LineWidth',2)
xlim([-1.5*delta,1.5*delta])
ylim([0,1])
hold on
plot([delta delta],[0 1],'k--',[-delta -delta],[0 1],'k--','LineWidth',1)
hold off
legend('Power_{H}','Power_{H} 2\times H','Power_{P}','Power_{P} 2\times P','Location','North')
xlabel('w')
ylabel('Power')
UtilStdizeFigureAbs(2);

% Save the figure if dosave = true
if nargin ~= 1
    dosave = false;
end
if dosave
    UtilSaveFigFile(1, 'Profher','ProfherPower_CPCS', [],'eps');
    UtilSaveFigFile(2, 'Profher','ProfherPower_Power', [],'eps');
end