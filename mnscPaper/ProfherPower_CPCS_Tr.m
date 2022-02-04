function ProfherPower_CPCS_Tr(dosave)
close all
basic = ProfherSetParameters();
[astar,~] = OneShotMaxExpectedNetGain(basic,true,true);
delta = basic.relevantdifference;
w = linspace(-1.5*delta,1.5*delta,delta/2);
CPCSmkt = OneShotCPCD(w,basic,astar(1),astar(2));
Powermkt = PowerPlotTwoSided(w,basic,astar(1),astar(2),'Frequentist','Frequentist');

basic = MixInputModifier(basic,{'horizon','Patent'});
[astar,~] = OneShotMaxExpectedNetGain(basic,true,true);
CPCSpat = OneShotCPCD(w,basic,astar(1),astar(2));
Powerpat = PowerPlotTwoSided(w,basic,astar(1),astar(2),'Frequentist','Frequentist');

CPCSactual = OneShotCPCD(w,basic,32,7.833);
Poweractual = PowerPlotTwoSided(w,basic,32,7.833,'Frequentist','Frequentist');

figure(1)
plot(w,CPCSpat,'b-',w,CPCSmkt,'k--',w,CPCSactual,'r-.','LineWidth',2)
xlim([-1.5*delta,1.5*delta])
ylim([0,1])
hold on
plot([delta delta],[0 1],'k:',[-delta -delta],[0 1],'k:','LineWidth',1)
hold off
legend('CPCS_{H}','CPCS_{P}','CPCS actual trial','Location','South')
xlabel('W')
ylabel('CPCS')
UtilStdizeFigureAbs(1);

figure(2)
plot(w,Powerpat,'b-',w,Powermkt,'k--',w,Poweractual,'r-.','LineWidth',2)
xlim([-1.5*delta,1.5*delta])
ylim([0,1])
hold on
plot([delta delta],[0 1],'k:',[-delta -delta],[0 1],'k:','LineWidth',1)
hold off
legend('Power_{H}','Power_{P}','Power actual trial','Location','North')
xlabel('W')
ylabel('Power')
UtilStdizeFigureAbs(2);

% Save the figure if dosave = true
if nargin ~= 1
    dosave = false;
end
if dosave
    UtilSaveFigFile(1, 'Profher','ProfherPower_CPCS_Tr', [],'eps');
    UtilSaveFigFile(2, 'Profher','ProfherPower_Power_Tr', [],'eps');
end