function ProfherWTPrecrate(dosave)
close all
basic = ProfherSetParameters();
cr = basic.cr;
basic = MixInputModifier(basic,{'cfix',0,'cr',0});
r0 = basic.r;
zeta = basic.zeta;
h = 0.3;
r = 1:h:50;

[~,ENGstar,~,~] = SensitivityAnalysisOneParameter(basic,'r',r);
ENGstarDiff = diff(ENGstar)./h;
figure(1)
semilogy(r(1:end-1),ENGstarDiff,'b-','LineWidth',2)
figure(2)
plot(r,ENGstar,'b-','LineWidth',2)

basic = MixInputModifier(basic,{'horizon','Patent'});
[~,ENGstar,~,~] = SensitivityAnalysisOneParameter(basic,'r',r);
ENGstarDiff = diff(ENGstar)./h;
figure(1)
hold on
plot(r(1:end-1),ENGstarDiff,'r--','LineWidth',2)
figure(2)
hold on
plot(r,ENGstar,'r--','LineWidth',2)

basic = MixInputModifier(basic,{'CTdiscrate',0});
[~,ENGstar,~,~] = SensitivityAnalysisOneParameter(basic,'r',r);
ENGstarDiff = diff(ENGstar)./h;
figure(1)
hold on
plot(r(1:end-1),ENGstarDiff,'r-.','LineWidth',2)
plot(r,cr.*ones(size(r)),'k:','LineWidth',2);
xlabel('Recruitment rate r (mo^{-1})')
ylabel('Marginal value (cost) of an additional unit of r')
legend('mv_{r,P} discounted','mv_{r,H} discounted','mv_{r,H} undiscounted','c_r')
hold off
UtilStdizeFigureAbs(1);
figure(2)
hold on
plot(r,ENGstar,'r-.','LineWidth',2)
xlabel('Recruitment rate r (mo^{-1})')
ylabel('Optimal expected Net Gain (£)')
legend('V^*_P discounted','V*_H discounted','V*_H undiscounted','Location','SouthEast')
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
