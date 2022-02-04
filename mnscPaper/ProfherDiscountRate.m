function ProfherDiscountRate(dosave)
close all
basic = ProfherSetParameters();

rho = 0.0005:0.0005:0.01;
rho_per_annum_perc = (exp(12*rho) - 1)*100;
rho_base = (exp(12*basic.CTdiscrate) - 1)*100;
rng(121);

for j = 1:length(rho)
    basic = MixInputModifier(basic,{'CTdiscrate',rho(j)});
    [Varstarmkt(:,j),ENGstarmkt(j)] = OneShotMaxExpectedNetGain(basic,true,true);
end
Tstarmkt = Varstarmkt(1,:);
rstarmkt = Varstarmkt(2,:);
Qstarmkt = Tstarmkt .* rstarmkt./2;


basic = MixInputModifier(basic,{'horizon','Patent'});
for j = 1:length(rho)
    basic = MixInputModifier(basic,{'CTdiscrate',rho(j)});
    [Varstarpat(:,j),ENGstarpat(j)] = OneShotMaxExpectedNetGain(basic,true,true);
end
Tstarpat = Varstarpat(1,:);
rstarpat = Varstarpat(2,:);
Qstarpat = Tstarpat .* rstarpat./2;

figure(1)
plot(rho_per_annum_perc,Qstarpat,'b-',rho_per_annum_perc,Qstarmkt,'k--','LineWidth',2)
hold on
xlabel('Annual discount rate (%)');
ylabel('Optimal sample size Q^*')
ax = axis;
ax(3) = 0;
axis(ax)
plot([rho_base rho_base],[ax(3) ax(4)],'k:','LineWidth',2)
legend('Q^*_{H}','Q^*_{P}','Base parameter','Location','NorthEast')
UtilStdizeFigureAbs(1);

figure(2)
plot(rho_per_annum_perc,Tstarpat,'b-',rho_per_annum_perc,Tstarmkt,'k--','LineWidth',2)
hold on
xlabel('Annual discount rate (%)');
ylabel('Optimal recruiment period duration T^* (mo)')
ax = axis;
ax(3) = 0;
axis(ax)
plot([rho_base rho_base],[ax(3) ax(4)],'k:','LineWidth',2)
legend('T^*_{H}','T^*_{P}','Base parameter','Location','NorthEast')
UtilStdizeFigureAbs(2);

figure(3)
plot(rho_per_annum_perc,rstarpat,'b-',rho_per_annum_perc,rstarmkt,'k--','LineWidth',2)
hold on
xlabel('Annual discount rate (%)');
ylabel('Optimal recruitment rate r^* (mo^{-1})')
ax = axis;
ax(3) = 0;
axis(ax)
plot([rho_base rho_base],[ax(3) ax(4)],'k:','LineWidth',2)
legend('r^*_{H}','r^*_{P}','Base parameter','Location','SouthEast')
UtilStdizeFigureAbs(3);

figure(4)
plot(rho_per_annum_perc,ENGstarpat./1e6,'b-',rho_per_annum_perc,ENGstarmkt./1e6,'k--','LineWidth',2)
hold on
ax = axis;
ax(3) = 0;
axis(ax)
plot([rho_base rho_base],[ax(3) ax(4)],'k:','LineWidth',2)
xlabel('Annual discount rate (%)');
ylabel('Optimal expected net gain V^* (millions of £)')
legend('V^*_{H}','V^*_{P}','Base parameter','Location','NorthEast')
UtilStdizeFigureAbs(4);


% Save the figures if dosave=true
if nargin < 1
    dosave = false;
end
if dosave
    UtilSaveFigFile(1, 'Profher', 'ProfherDiscrate_Tr_Q',[],'eps');
    UtilSaveFigFile(2, 'Profher', 'ProfherDiscrate_Tr_T',[],'eps');
    UtilSaveFigFile(3, 'Profher', 'ProfherDiscrate_Tr_r',[],'eps');
    UtilSaveFigFile(4, 'Profher', 'ProfherDiscrate_Tr_ENG',[],'eps');
end

