function ProfherOnlineOffline(dosave)
close all

basic = ProfherSetParameters();
mult = 1/1;
basic = MixInputModifier(basic,{'Population',basic.Population*mult,'H',basic.H*mult});
online = [0 1];
delta = basic.relevantdifference;
mu0 = -2*delta:(delta/10):2*delta;
rng(122);

basic = MixInputModifier(basic,{'online',online(1)});
for j = 1:length(mu0)
    basic = MixInputModifier(basic,{'mu0',mu0(j)});
    [Varstarmkt(:,j),ENGstarmkt1(j)] = OneShotMaxExpectedNetGain(basic,true,'Tr',300);
end
Tstarmkt1 = Varstarmkt(1,:);
rstarmkt1 = Varstarmkt(2,:);
Qstarmkt1 = Tstarmkt1 .* rstarmkt1./2;


basic = MixInputModifier(basic,{'horizon','Patent'});
for j = 1:length(mu0)
    basic = MixInputModifier(basic,{'mu0',mu0(j)});
    [Varstarpat(:,j),ENGstarpat1(j)] = OneShotMaxExpectedNetGain(basic,true,'Tr',300);
end
Tstarpat1 = Varstarpat(1,:);
rstarpat1 = Varstarpat(2,:);
Qstarpat1 = Tstarpat1 .* rstarpat1./2;

basic = MixInputModifier(basic,{'horizon','MktExcl','online',online(2)});
for j = 1:length(mu0)
    basic = MixInputModifier(basic,{'mu0',mu0(j)});
    [Varstarmkt(:,j),ENGstarmkt2(j)] = OneShotMaxExpectedNetGain(basic,true,'Tr',300);
end
Tstarmkt2 = Varstarmkt(1,:);
rstarmkt2 = Varstarmkt(2,:);
Qstarmkt2 = Tstarmkt2 .* rstarmkt2./2;


basic = MixInputModifier(basic,{'horizon','Patent'});
for j = 1:length(mu0)
    basic = MixInputModifier(basic,{'mu0',mu0(j)});
    [Varstarpat(:,j),ENGstarpat2(j)] = OneShotMaxExpectedNetGain(basic,true,'Tr',300);
end
Tstarpat2 = Varstarpat(1,:);
rstarpat2 = Varstarpat(2,:);
Qstarpat2 = Tstarpat2 .* rstarpat2./2;

figure(1)
plot(mu0,Qstarpat1,'b-',mu0,Qstarpat2,'b:',mu0,Qstarmkt1,'k--',mu0,Qstarmkt2,'k-.','LineWidth',2)
xlabel('\mu_0');ylabel('Optimal sample size Q^*')
lgd = legend('Q^*_{H} offline','Q^*_{H} online','Q^*_{P} offline','Q^*_{P} online','Location','South');
ax = axis;
ax(3) = 0;
axis(ax)
UtilStdizeFigureAbs(1);

figure(2)
plot(mu0,Tstarpat1,'b-',mu0,Tstarpat2,'b:',mu0,Tstarmkt1,'k--',mu0,Tstarmkt2,'k-.','LineWidth',2)
xlabel('\mu_0');ylabel('Optimal recruitment duration T^*')
lgd = legend('T^*_{H} offline','T^*_{H} online','T^*_{P} offline','T^*_{P} online','Location','South');
ax = axis;
ax(3) = 0;
axis(ax)
UtilStdizeFigureAbs(2);

figure(3)
plot(mu0,rstarpat1,'b-',mu0,rstarpat2,'b:',mu0,rstarmkt1,'k--',mu0,rstarmkt2,'k-.','LineWidth',2)
xlabel('\mu_0');ylabel('Optimal recruitment rate r^*')
lgd = legend('r^*_{H} offline','r^*_{H} online','r^*_{P} offline','r^*_{P} online','Location','South');
ax = axis;
ax(3) = 0;
axis(ax)
UtilStdizeFigureAbs(3);

figure(4)
plot(mu0,ENGstarpat1/1e6,'b-',mu0,ENGstarpat2/1e6,'b:',mu0,ENGstarmkt1/1e6,'k--',mu0,ENGstarmkt2/1e6,'k-.','LineWidth',2)
xlabel('\mu_0');ylabel('Optimal expected net gain V^* (millions of £)')
lgd = legend('V^*_{H} offline','V^*_{H} online','V^*_{P} offline','V^*_{P} online','Location','SouthEast');
ax = axis;
ax(3) = 0;
axis(ax)
UtilStdizeFigureAbs(4);


% Save the figures if dosave=true
if nargin < 1
    dosave = false;
end
if dosave
    UtilSaveFigFile(1, 'Profher', 'ProfherOnOff_Q',[],'eps');
    UtilSaveFigFile(2, 'Profher', 'ProfherOnOff_T',[],'eps');
    UtilSaveFigFile(3, 'Profher', 'ProfherOnOff_r',[],'eps');
    UtilSaveFigFile(4, 'Profher', 'ProfherOnOff_ENG',[],'eps');
end



