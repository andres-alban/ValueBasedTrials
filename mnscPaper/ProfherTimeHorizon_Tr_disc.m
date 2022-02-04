function ProfherTimeHorizon_Tr_disc(dosave)

close all
basic = ProfherSetParameters();
H_base = basic.H;

H = [12,15,18,20:20:400];
P = basic.zeta.*H;
rng(122);

for j = 1:length(P)
    basic = MixInputModifier(basic,{'Population',P(j)});
    [Varstarmkt(:,j),ENGstarmkt(j)] = OneShotMaxExpectedNetGain(basic,true,true,200);
end
Tstarmkt = Varstarmkt(1,:);
rstarmkt = Varstarmkt(2,:);
Qstarmkt = Tstarmkt .* rstarmkt./2;
[Tasympt,rasympt,ENGasympt] = OneShotTrasymptoticapprox(basic,P);
Qasympt = Tasympt .* rasympt./2;


basic = MixInputModifier(basic,{'horizon','Patent'});
for j = 1:length(H)
    basic = MixInputModifier(basic,{'H',H(j)});
    [Varstarpat(:,j),ENGstarpat(j)] = OneShotMaxExpectedNetGain(basic,true,true,200);
end
Tstarpat = Varstarpat(1,:);
rstarpat = Varstarpat(2,:);
Qstarpat = Tstarpat .* rstarpat./2;

figure(1)
plot(H/12,Qstarpat,'b-',H/12,Qstarmkt,'k--',H/12,Qasympt,'r-.','LineWidth',2)
hold on
plot([H_base/12 H_base/12],[0 Qasympt(1)],'k:','LineWidth',2)
xlabel('Time horizon H (year)');
ylabel('Optimal sample size Q^*')
legend('Q^*_{H}','Q^*_{P}','Asymptote','Base parameter','Location','SouthEast')
UtilStdizeFigureAbs(1);

figure(2)
plot(H/12,Tstarpat,'b-',H/12,Tstarmkt,'k--',H/12,Tasympt,'r-.','LineWidth',2)
hold on
plot([H_base/12 H_base/12],[0 Tstarmkt(1)],'k:','LineWidth',2)
xlabel('Time horizon H (year)');
ylabel('Optimal recruiment period duration T^*(mo)')
legend('T^*_{H}','T^*_{P}','Asymptote','Base parameter','Location','SouthEast')
UtilStdizeFigureAbs(2);

figure(3)
plot(H/12,rstarpat,'b-',H/12,rstarmkt,'k--',H/12,rasympt,'r-.','LineWidth',2)
hold on
plot([H_base/12 H_base/12],[0 rasympt(1)],'k:','LineWidth',2)
xlabel('Time horizon H (year)');
ylabel('Optimal recruitment rate r^* (mo^{-1})')
legend('r^*_{H}','r^*_{P}','Asymptote','Base parameter','Location','SouthEast')
UtilStdizeFigureAbs(3);

figure(4)
plot(H/12,ENGstarpat/10e6,'b-',H/12,ENGstarmkt/10e6,'k--',H/12,ENGasympt/10e6,'r-.','LineWidth',2)
hold on
plot([H_base/12 H_base/12],[0 ENGasympt(1)/10e6],'k:','LineWidth',2)
xlabel('Time horizon H (year)');
ylabel('Optimal expected net gain V^* (millions of £)')
legend('V^*_{H}','V^*_{P}','Asymptote','Base parameter','Location','SouthEast')
UtilStdizeFigureAbs(4);


% Save the figures if dosave=true
if nargin < 1
    dosave = false;
end
if dosave
    UtilSaveFigFile(1, 'Profher', 'ProfherTimeHorizon_Tr_Q_new',[],'eps');
    UtilSaveFigFile(2, 'Profher', 'ProfherTimeHorizon_Tr_T_new',[],'eps');
    UtilSaveFigFile(3, 'Profher', 'ProfherTimeHorizon_Tr_r_new',[],'eps');
    UtilSaveFigFile(4, 'Profher', 'ProfherTimeHorizon_Tr_ENG_new',[],'eps');
end








