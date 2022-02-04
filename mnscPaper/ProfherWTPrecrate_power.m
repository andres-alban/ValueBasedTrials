function ProfherWTPrecrate_power(dosave)
close all
basic = ProfherSetParameters();
cr = basic.cr;
basic = MixInputModifier(basic,{'cfix',0,'cr',0});
r0 = basic.r;
zeta = basic.zeta;
h = 0.3;
r = 1:h:50;

basic = MixInputModifier(basic,{'horizon','Patent'});
[~,ENGstar,~,~] = SensitivityAnalysisOneParameter(basic,'r',r);
ENGstarDiff = diff(ENGstar./1e6)./h;
figure(1)
semilogy(r(1:end-1),ENGstarDiff,'b-','LineWidth',2)
figure(2)
plot(r,ENGstar./1e6,'b-','LineWidth',2)

basic = MixInputModifier(basic,{'horizon','MktExcl'});
[~,ENGstar,~,~] = SensitivityAnalysisOneParameter(basic,'r',r);
ENGstarDiff = diff(ENGstar./1e6)./h;
figure(1)
hold on
semilogy(r(1:end-1),ENGstarDiff,'k--','LineWidth',2)
figure(2)
hold on
plot(r,ENGstar./1e6,'k--','LineWidth',2)

figure(1)
hold on
% functional form of the cost per recruitment rate
% Adjust marginal cost according to the output of ProfherSetup_cost
marginal_cost = 3.06*1.37e4*0.0559*r.^2.06;
plot(r,marginal_cost./1e6,'r:','LineWidth',2);
xlabel('Recruitment rate r (mo^{-1})')
ylabel({'Marginal value and cost of an additional unit of r','(millions of £)'})
legend('mv_{r,H}','mv_{r,P}','mc_r')
hold off
UtilStdizeFigureAbs(1);
figure(2)
xlabel('Recruitment rate r (mo^{-1})')
ylabel('Optimal expected Net Gain (£)')
legend('V^*_H','V*_P','Location','SouthEast')
hold off
UtilStdizeFigureAbs(2);


% Save the figures if dosave=true
if nargin < 1
    dosave = false;
end
if dosave
    UtilSaveFigFile(1, 'Profher', 'ProfherWTPrecrate',[],'eps');
    UtilSaveFigFile(2, 'Profher', 'Profherr_B',[],'eps');
end
