function ProfherCostAddSite(horizon,online,dosave)

close all
% If not enough inputs the default is online learning and patent protection
if nargin < 2
    online = 1;
end
if nargin < 1
    horizon = 'Patent';
end

if strcmp(horizon,'MktExcl')
    subhorizon = 'mkt';
    titlehorizon = 'Market Exclusivity';
elseif strcmp(horizon,'Patent')
    subhorizon = 'pat';
    titlehorizon = 'Patent Protection';
else
    error('horizon wrongly specified. Choose MktExcl or Patent')
end

if online == 1
    subonline = 'on';
    titleonline = 'Online learning';
elseif online == 0
    subonline = 'off';
    titleonline = 'Offline learning';
else
    error('online wrongly specified. Choose 1 for online or 0 for offline')
end
subscript = ['{',subonline,',',subhorizon,'}'];
% Define the model with the appropriate parameters for Profher
basic = ProfherSetParameters();
basic = MixInputModifier(basic,{'horizon',horizon,'online',online});

zeta = basic.zeta;
r = 20:1:1002;
mu0 = -1500:1500:1500;
plotstyle = {'r:','b-','k--'};
leg = cell(size(mu0));
for i = 1:length(mu0)
    basic = MixInputModifier(basic,{'mu0',mu0(i)});
    [~,ENGstar,~,p] = SensitivityAnalysisOneParameter(basic,'r',r);
    ENGstarDiff = diff(ENGstar);
    figure(1)
    hold on
    plot(p(1:end-1),ENGstarDiff,plotstyle{i},'LineWidth',3)
    leg{i} = sprintf('\\mu_0 = %d',mu0(i));
end
hold off
xlabel('Recruitment rate')
ylabel('Willingness to pay for an additional unit of recruitment')
title([titleonline,' - ',titlehorizon])
legend(leg)
UtilStdizeFigureAbs(1);

% Save the figures if dosave=true
if nargin <3
    dosave = false;
end
if dosave
    UtilSaveFigFile(1, 'Profher', 'WTPaddSite',[subonline,horizon],'eps');
end

