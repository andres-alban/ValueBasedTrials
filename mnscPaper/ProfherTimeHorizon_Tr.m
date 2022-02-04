function ProfherTimeHorizon_Tr(dosave)
% Outdated
close all
basic = ProfherSetParameters();
basic = MixInputModifier(basic,{'cfix',basic.cfix});

H = 100:20:300;
P = basic.zeta.*H;
rng(121);

for j = 1:length(P)
    basic = MixInputModifier(basic,{'Population',P(j)});
    [Varstarmkt(:,j),ENGstarmkt(j)] = OneShotMaxExpectedNetGain(basic,true,true);
end
Tstarmkt = Varstarmkt(1,:);
rstarmkt = Varstarmkt(2,:);
Qstarmkt = Tstarmkt .* rstarmkt./2;
[Tasympt,rasympt,ENGasympt] = OneShotTrasymptoticapprox(basic,P);
Qasympt = Tasympt .* rasympt./2;


basic = MixInputModifier(basic,{'horizon','Patent'});
for j = 1:length(H)
    basic = MixInputModifier(basic,{'H',H(j)});
    [Varstarpat(:,j),ENGstarpat(j)] = OneShotMaxExpectedNetGain(basic,true,true);
end
Tstarpat = Varstarpat(1,:);
rstarpat = Varstarpat(2,:);
Qstarpat = Tstarpat .* rstarpat./2;


basic = MixInputModifier(basic,{'CTdiscrate',0});
for j = 1:length(H)
    basic = MixInputModifier(basic,{'H',H(j)});
    [Varstarpatund(:,j),ENGstarpatund(j)] = OneShotMaxExpectedNetGain(basic,true,true);
end
Tstarpatund = Varstarpatund(1,:);
rstarpatund = Varstarpatund(2,:);
Qstarpatund = Tstarpatund .* rstarpatund./2;
[Tasymptpatund,rasymptpatund,ENGasymptpatund] = OneShotTrasymptoticapprox(basic,H);
Qasymptpatund = Tasymptpatund .* rasymptpatund./2;


figure(1)
plot(H,Qstarpatund,'b-',H,Qstarpat,'b-.',H,Qstarmkt,'k--',H,Qasymptpatund,'rx-',H,Qasympt,'r:','LineWidth',2)
xlabel('Time horizon H (mo)');
ylabel('Optimal number of pairwise allocations Q^*')
ax = axis;
ax(4) = ax(4) + 80;
axis(ax)
legend('Q^*_{H} undiscounted','Q^*_{H} discounted','Q^*_{P} discounted','Asymptotic approx. Q^*_{H} undiscounted','Asymptotic approx. discounted','Location','NorthWest')
UtilStdizeFigureAbs(1);

figure(2)
plot(H,Tstarpatund,'b-',H,Tstarpat,'b-.',H,Tstarmkt,'k--',H,Tasymptpatund,'rx-',H,Tasympt,'r:','LineWidth',2)
xlabel('Time horizon H (mo)');
ylabel('Optimal recruiment period duration T^*')
ax = axis;
ax(4) = ax(4) + 2;
axis(ax)
legend('T^*_{H} undiscounted','T^*_{H} discounted','T^*_{P} discounted','Asymptotic approx. T^*_{H} undiscounted','Asymptotic approx. discounted','Location','NorthWest')
UtilStdizeFigureAbs(2);

figure(3)
plot(H,rstarpatund,'b-',H,rstarpat,'b-.',H,rstarmkt,'k--',H,rasymptpatund,'rx-',H,rasympt,'r:','LineWidth',2)
xlabel('Time horizon H (mo)');
ylabel('Optimal recruitment rate r^*')
ax = axis;
ax(4) = ax(4) + 50;
axis(ax)
legend('r^*_{H} undiscounted','r^*_{H} discounted','r^*_{P} discounted','Asymptotic approx. r^*_{H} undiscounted','Asymptotic approx. discounted','Location','NorthWest')
UtilStdizeFigureAbs(3);

figure(4)
plot(H,ENGstarpatund,'b-',H,ENGstarpat,'b-.',H,ENGstarmkt,'k--',H,ENGasymptpatund,'rx-',H,ENGasympt,'r:','LineWidth',2)
xlabel('Time horizon H (mo)');
ylabel('Optimal expected net gain V^*')
legend('V^*_{H} undiscounted','V^*_{H} discounted','V^*_{P} discounted','Asymptotic approx. V^*_{H} undiscounted','Asymptotic approx. discounted','Location','NorthWest')
UtilStdizeFigureAbs(4);


% Save the figures if dosave=true
if nargin < 1
    dosave = false;
end
if dosave
    UtilSaveFigFile(1, 'Profher', 'ProfherTimeHorizon_Tr_Q',[],'eps');
    UtilSaveFigFile(2, 'Profher', 'ProfherTimeHorizon_Tr_T',[],'eps');
    UtilSaveFigFile(3, 'Profher', 'ProfherTimeHorizon_Tr_r',[],'eps');
    UtilSaveFigFile(4, 'Profher', 'ProfherTimeHorizon_Tr_ENG',[],'eps');
end








