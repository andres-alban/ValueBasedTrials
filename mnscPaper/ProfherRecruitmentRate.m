function ProfherRecruitmentRate(dosave)
close all
basic = ProfherSetParameters();

r = 1:1:25;
rng(121);

[VarstarTrue,ENGstarmktTrue] = OneShotMaxExpectedNetGain(basic,true,'Tr',100,false);
TstarmktTrue = VarstarTrue(1);
rstarmktTrue = VarstarTrue(2);
QstarmktTrue = TstarmktTrue .* rstarmktTrue./2;
for j = 1:length(r)
    basic = MixInputModifier(basic,{'r',r(j)});
    [Varstarmkt(j),ENGstarmkt(j)] = OneShotMaxExpectedNetGain(basic,true,'T',20,false);
    ENGmkt(j) = OneShotExpectedNetGain(basic,Varstarmkt(j),r(j));
end
Tstarmkt = Varstarmkt;
% rstarmkt = Varstarmkt(2,:);
Qstarmkt = Tstarmkt .* r./2;


basic = MixInputModifier(basic,{'horizon','Patent'});
[VarstarTrue,ENGstarpatTrue] = OneShotMaxExpectedNetGain(basic,true,'Tr',20,false);
TstarpatTrue = VarstarTrue(1);
rstarpatTrue = VarstarTrue(2);
QstarpatTrue = TstarpatTrue .* rstarpatTrue./2;
for j = 1:length(r)
    basic = MixInputModifier(basic,{'r',r(j)});
    [Varstarpat(j),ENGstarpat(j)] = OneShotMaxExpectedNetGain(basic,true,'T',100,false);
    ENGpat(j) = OneShotExpectedNetGain(basic,Varstarpat(j),r(j));
end
Tstarpat = Varstarpat;
% rstarpat = Varstarpat(2,:);
Qstarpat = Tstarpat .* r./2;

figure(1)
plot(r,Qstarpat,'b-',r,ones(size(Qstarpat)).*QstarpatTrue,'b:',r,Qstarmkt,'k--',r,ones(size(Qstarmkt)).*QstarmktTrue,'k-.','LineWidth',2)
hold on
xlabel('Recruitment rate (mo^{-1})');
ylabel('Sample size Q')
ax = axis;
ax(3) = 0;
axis(ax)
% plot([rstarmktTrue rstarmktTrue],[ax(3) ax(4)],'k:','LineWidth',2)
legend('Q^*_{H}','Q^*_{H}','Q^*_{P}','Q^*_{P}','Location','NorthEast')
UtilStdizeFigureAbs(1);



% figure(2)
% plot(r,Tstarpat,'b-',r,Tstarmkt,'k--','LineWidth',2)
% hold on
% xlabel('Recruitment rate (mo^{-1})');
% ylabel('Optimal recruiment period duration T^* (mo)')
% ax = axis;
% ax(3) = 0;
% axis(ax)
% plot([rho_base rho_base],[ax(3) ax(4)],'k:','LineWidth',2)
% legend('T^*_{H}','T^*_{P}','Base parameter','Location','NorthEast')
% UtilStdizeFigureAbs(2);

% figure(3)
% plot(rho_per_annum_perc,rstarpat,'b-',rho_per_annum_perc,rstarmkt,'k--','LineWidth',2)
% hold on
% xlabel('Annual discount rate (%)');
% ylabel('Optimal recruitment rate r^* (mo^{-1})')
% ax = axis;
% ax(3) = 0;
% axis(ax)
% plot([rho_base rho_base],[ax(3) ax(4)],'k:','LineWidth',2)
% legend('r^*_{H}','r^*_{P}','Base parameter','Location','SouthEast')
% UtilStdizeFigureAbs(3);

figure(4)
plot(r,ENGpat./1e6,'b-',r,ones(size(ENGpat)).*ENGstarpatTrue./1e6,'b:',r,ENGmkt./1e6,'k--',r,ones(size(ENGpat)).*ENGstarmktTrue./1e6,'k-.','LineWidth',2)
hold on
ax = axis;
ax(3) = 0;
axis(ax)
% plot([rstarmktTrue rstarmktTrue],[ax(3) ax(4)],'k:','LineWidth',2)
xlabel('Recruitment rate (mo^{-1})');
ylabel('Expected net gain V (millions of £)')
legend('V_{H}','V^*_{H}','V_{P}','V^*_{P}','Location','NorthEast')
UtilStdizeFigureAbs(4);

figure(5)
plot(r,(1 - ENGpat./(ones(size(ENGpat)).*ENGstarpatTrue)).*100,'b-',r,(1 - ENGmkt./(ones(size(ENGpat)).*ENGstarmktTrue)).*100,'k--','LineWidth',2)
hold on
ax = axis;
ax(3) = 0;
ax(4) = 40;
axis(ax)
xlabel('Recruitment rate (mo^{-1})');
ylabel('Percentage loss in expected net gain')
legend('Fixed horizon','Fixed pool','Location','NorthEast')
UtilStdizeFigureAbs(5);


% Save the figures if dosave=true
if nargin < 1
    dosave = false;
end
if dosave
    UtilSaveFigFile(1, 'Profher', 'ProfherRecrate_Q',[],'eps');
    UtilSaveFigFile(4, 'Profher', 'ProfherRecrate_ENG',[],'eps');
    UtilSaveFigFile(5, 'Profher', 'ProfherRecrate_ENG_relative',[],'eps');
end

