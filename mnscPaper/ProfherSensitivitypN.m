function ProfherSensitivitypN(horizon,online,dosave)
% This function evaluates the effect of pN on Tstar nad ENGstar for
% Profher. It plots some of the figures in section 7 of the paper.
% 
%  Revision: AA 26-11-2017

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
elseif strcmp(horizon,'Patent')
    subhorizon = 'pat';
else
    error('horizon wrongly specified. Choose MktExcl or Patent')
end

if online == 1
    subonline = 'on';
elseif online == 0
    subonline = 'off';
else
    error('online wrongly specified. Choose 1 for online or 0 for offline')
end
subscript = ['{',subonline,',',subhorizon,'}'];
% Define the model with the appropriate parameters for Profher
basic = ProfherSetParameters();
basic = MixInputModifier(basic,{'horizon',horizon,'online',online});

% Plot Tstar and ENGstar as function of pN
pN = 0.01:0.01:0.99;
mu0 = -1500:1500:1500;
plotstyle = {'r:','b-','k--'};
leg = cell(size(mu0));
for i = 1:length(mu0)
    basic = MixInputModifier(basic,{'mu0',mu0(i)});
    [Tstar,ENGstar,PCDstar,p] = SensitivityAnalysisOneParameter(basic,'pN',pN);
    Qstar = basic.r .* Tstar;
    figure(1)
    hold on
    plot(p,Tstar,plotstyle{i},'LineWidth',3)
    figure(2)
    hold on
    plot(p,ENGstar,plotstyle{i},'LineWidth',3)
    leg{i} = sprintf('\\mu_0 = %.2f',mu0(i));
    figure(3)
    hold on
    plot(p,Qstar,plotstyle{i},'LineWidth',3)
    figure(4)
    hold on
    plot(p,PCDstar,plotstyle{i},'LineWidth',3)
    leg{i} = sprintf('\\mu_0 = %d',mu0(i));
end
figure(1)
xlabel('Proportion currently in new treatment (p_N)');ylabel(['Optimal trial length T^*_',subscript])
legend(leg,'Location','North')
hold off
UtilStdizeFigureAbs(1);
figure(2)
xlabel('Proportion currently in new treatment (p_N)');ylabel(['Maximum expected net gain B^*_',subscript])
legend(leg,'Location','North')
hold off
UtilStdizeFigureAbs(2);
figure(3)
xlabel('Proportion currently in new treatment (p_N)');ylabel(['Optimal number of pairwise allocations Q^*_',subscript])
legend(leg,'Location','North')
hold off
UtilStdizeFigureAbs(3);
figure(4)
xlabel('Proportion currently in new treatment (p_N)');ylabel(['Probability of correct decision PCD_',subscript])
legend(leg,'Location','North')
hold off
UtilStdizeFigureAbs(4);

%%%%%%%%%%%%%%%%%%%%%
% Power curves
%%%%%%%%%%%%%%%%%%%%%
basic = MixInputModifier(basic,{'mu0',0});
pN = [0.3 0.5 0.7];
delta = basic.relevantdifference;
w = linspace(-1.5*delta,1.5*delta,delta);
plotstyle = {'r:','b-','k--'};
leg = cell(size(pN));
figure(5)
hold on
for i = 1:length(pN)
    basic = MixInputModifier(basic,{'pN',pN(i)});
    [Tstar,~] = OneShotMaxExpectedNetGain(basic);
    power = OneShotCPCD(w,basic,Tstar);
    plot(w,power,plotstyle{i},'LineWidth',3)
    leg{i} = sprintf('p_N = %0.1f',pN(i));
end
ylim([0 1])
legend(leg)
xlabel('True expected INMB, W')
ylabel('Power')
hold off
UtilStdizeFigureAbs(5);
basic = MixInputModifier(basic,{'pN',0.39});

% Save the figures if dosave is true
if nargin < 3
    dosave = false;
end

if dosave
    UtilSaveFigFile(1, 'Profher', 'OptTrialLength_vs_pN',[subonline,horizon],'eps');
    UtilSaveFigFile(2, 'Profher', 'MaxENG_vs_pN',[subonline,horizon],'eps');
    UtilSaveFigFile(3, 'Profher', 'OptNumPairs_vs_pN',[subonline,horizon],'eps');
    UtilSaveFigFile(4, 'Profher', 'PCD_vs_pN',[subonline,horizon],'eps');
    UtilSaveFigFile(5, 'Profher', 'CPCD_vs_pN',[subonline,horizon],'eps');

end

