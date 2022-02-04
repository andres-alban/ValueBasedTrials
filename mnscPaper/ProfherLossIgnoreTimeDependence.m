function ProfherLossIgnoreTimeDependence(dosave)
close all
basic = ProfherSetParameters();
basicpat = MixInputModifier(basic,{'horizon','Patent'});
basic2 = basic;
basic2 = MixInputModifier(basic2,{'CTdiscrate',0,'cr',0});

r = 1:1:25;
rng(121);

[VarstarTrue,ENGstarmktTrue] = OneShotMaxExpectedNetGain(basic,true,'Tr',20,false);
TstarmktTrue = VarstarTrue(1);
rstarmktTrue = VarstarTrue(2);
QstarmktTrue = TstarmktTrue .* rstarmktTrue./2;

[VarstarTrue,ENGstarpatTrue] = OneShotMaxExpectedNetGain(basicpat,true,'Tr',20,false);
TstarpatTrue = VarstarTrue(1);
rstarpatTrue = VarstarTrue(2);
QstarpatTrue = TstarpatTrue .* rstarpatTrue./2;

[Varstar,ENGstar] = OneShotMaxExpectedNetGain(basic2,true,'T',20,false);
Qstar = Varstar * basic2.r / 2;

fprintf('Q^*_P = %.2f\n',QstarmktTrue)
fprintf('Q^*_H = %.2f\n',QstarpatTrue)
fprintf('Q^*_Ignore = %.2f\n',Qstar)

fprintf('ENG^*_P = %.2f Mill.\n',ENGstarmktTrue/1e6)
fprintf('ENG^*_H = %.2f Mill.\n',ENGstarpatTrue/1e6)
fprintf('ENG^*_ignore = %.2f Mill.\n',ENGstar/1e6)

ENGmkt = zeros(length(r),1);
ENGpat = zeros(length(r),1);
for j = 1:length(r)
    ENGmkt(j) = OneShotExpectedNetGain(basic,2*Qstar/r(j),r(j));
    ENGpat(j) = OneShotExpectedNetGain(basicpat,2*Qstar/r(j),r(j));
end


figure(1)
plot(r,ENGpat./1e6,'b-',r,ones(size(ENGpat)).*ENGstarpatTrue./1e6,'b:',r,ENGmkt./1e6,'k--',r,ones(size(ENGpat)).*ENGstarmktTrue./1e6,'k-.','LineWidth',2)
hold on
ax = axis;
ax(3) = 0;
axis(ax)
xlabel('Recruitment rate (mo^{-1})');
ylabel('Expected net gain V (millions of £)')
legend('V_{H}','V^*_{H}','V_{P}','V^*_{P}','Location','SouthEast')
UtilStdizeFigureAbs(1);

figure(2)
plot(r,(1 - ENGpat./(ones(size(ENGpat)).*ENGstarpatTrue)).*100,'b-',r,(1 - ENGmkt./(ones(size(ENGpat)).*ENGstarmktTrue)).*100,'k--','LineWidth',2)
hold on
ax = axis;
ax(3) = 0;
ax(4) = 40;
axis(ax)
xlabel('Recruitment rate (mo^{-1})');
ylabel('Percentage loss in expected net gain')
legend('Fixed horizon','Fixed pool','Location','NorthEast')
UtilStdizeFigureAbs(2);


% Save the figures if dosave=true
if nargin < 1
    dosave = false;
end
if dosave
    UtilSaveFigFile(1, 'Profher', 'ProfherIgnoreTime_ENG',[],'eps');
    UtilSaveFigFile(2, 'Profher', 'ProfherIgnoreTime_ENG_relative',[],'eps');
end

