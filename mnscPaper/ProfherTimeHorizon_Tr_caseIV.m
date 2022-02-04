function ProfherTimeHorizon_Tr_caseIV(dosave)

close all
basic = ProfherSetParameters();
basic = MixInputModifier(basic,{'cfix',basic.cfix});

H = [12,15,18,20:20:400];
P = basic.zeta.*H;
rng(121);

% for j = 1:length(P)
%     basic = MixInputModifier(basic,{'Population',P(j)});
%     [Varstarmkt(:,j),ENGstarmkt(j)] = OneShotMaxExpectedNetGain(basic,true,true,200);
% end
% Tstarmkt = Varstarmkt(1,:);
% rstarmkt = Varstarmkt(2,:);
% Qstarmkt = Tstarmkt .* rstarmkt./2;
% [Tasympt,rasympt,ENGasympt] = OneShotTrasymptoticapprox(basic,P);
% Qasympt = Tasympt .* rasympt./2;


basic = MixInputModifier(basic,{'horizon','Patent'});
% for j = 1:length(H)
%     basic = MixInputModifier(basic,{'H',H(j)});
%     [Varstarpat(:,j),ENGstarpat(j)] = OneShotMaxExpectedNetGain(basic,true,true,200);
% end
% Tstarpat = Varstarpat(1,:);
% rstarpat = Varstarpat(2,:);
% Qstarpat = Tstarpat .* rstarpat./2;


basic = MixInputModifier(basic,{'CTdiscrate',0});
for j = 1:length(H)
    basic = MixInputModifier(basic,{'H',H(j)});
    [Varstarpatund(:,j),ENGstarpatund(j)] = OneShotMaxExpectedNetGain(basic,true,true,200);
end
Tstarpatund = Varstarpatund(1,:);
rstarpatund = Varstarpatund(2,:);
Qstarpatund = Tstarpatund .* rstarpatund./2;
[Tasymptpatund,rasymptpatund,ENGasymptpatund] = OneShotTrasymptoticapprox(basic,H);
Qasymptpatund = Tasymptpatund .* rasymptpatund./2;


figure(1)
plot(H./12,Qstarpatund,'b-',H./12,Qasymptpatund,'r:','LineWidth',2)
xlabel('Time horizon H (year)');
ylabel('Optimal sample size Q^*')
ax = axis;
ax(4) = ax(4) + 80;
axis(ax)
legend('Q^*_{H} undiscounted','Asymptotic approx. Q^*_{H} undiscounted','Location','NorthWest')
UtilStdizeFigureAbs(1);

figure(2)
plot(H./12,Tstarpatund,'b-',H./12,Tasymptpatund,'r:','LineWidth',2)
xlabel('Time horizon H (year)');
ylabel('Optimal recruiment period duration T^* (mo)')
ax = axis;
ax(4) = ax(4) + 2;
axis(ax)
legend('T^*_{H} undiscounted','Asymptotic approx. T^*_{H} undiscounted','Location','NorthWest')
UtilStdizeFigureAbs(2);

figure(3)
plot(H./12,rstarpatund,'b-',H./12,rasymptpatund,'r:','LineWidth',2)
xlabel('Time horizon H (year)');
ylabel('Optimal recruitment rate r^* (mo^{-1})')
ax = axis;
ax(4) = ax(4) + 50;
axis(ax)
legend('r^*_{H} undiscounted','Asymptotic approx. r^*_{H} undiscounted','Location','NorthWest')
UtilStdizeFigureAbs(3);

figure(4)
plot(H./12,ENGstarpatund/1e6,'b-',H./12,ENGasymptpatund./1e6,'r:','LineWidth',2)
xlabel('Time horizon H (year)');
ylabel('Optimal expected net gain V^* (millions of £)')
legend('V^*_{H} undiscounted','Asymptotic approx. V^*_{H} undiscounted','Location','NorthWest')
UtilStdizeFigureAbs(4);


% Save the figures if dosave=true
if nargin < 1
    dosave = false;
end
if dosave
    UtilSaveFigFile(1, 'Profher', 'ProfherTimeHorizon_Tr_Q_undisc',[],'eps');
    UtilSaveFigFile(2, 'Profher', 'ProfherTimeHorizon_Tr_T_undisc',[],'eps');
    UtilSaveFigFile(3, 'Profher', 'ProfherTimeHorizon_Tr_r_undisc',[],'eps');
    UtilSaveFigFile(4, 'Profher', 'ProfherTimeHorizon_Tr_ENG_undisc',[],'eps');
end








