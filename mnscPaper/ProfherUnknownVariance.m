function ProfherUnknownVariance(horizon,online,dosave)
% This function compares the cases of known and unknown variance. 
% It plots some of the figues in section 7 of the paper.
%
% Revision: AA 19-12-2017

% If not enough inputs the default is online learning and patent protection

close all
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
basic = MixInputModifier(basic,{'online',online,'horizon',horizon});
basic = MixInputModifier(basic,{'UnkVariance',true});

% Plot Tstar and ENGstar as function of pN
xi0 = 2:0.1:20;
mu0 = -1500:1500:1500;
plotstyle = {'r:','b-','k--'};
leg = cell(size(mu0));
for i = 1:length(mu0)
    basic = MixInputModifier(basic,{'mu0',mu0(i)});
    Tstar = zeros(size(xi0));
    ENGstar = Tstar;
    PCDstar = Tstar;
    for j = 1:length(xi0)  % Loop over all the parameter values and find the maximum and maximizer of ENG
        basictemp = MixInputModifier(basic,{'UnkVarianceShape',xi0(j)});
        [Tstar(j),ENGstar(j)] = OneShotMaxExpectedNetGain(basictemp);
        PCDstar(j) = OneShotPCD(basictemp,Tstar(i));
    end
    Qstar = basic.r .* Tstar;
    figure(1)
    hold on
    plot(xi0,Tstar,plotstyle{i},'LineWidth',3)
    figure(2)
    hold on
    plot(xi0,ENGstar,plotstyle{i},'LineWidth',3)
    leg{i} = sprintf('\\mu_0 = %.2f',mu0(i));
    figure(3)
    hold on
    plot(xi0,Qstar,plotstyle{i},'LineWidth',3)
    figure(4)
    hold on
    plot(xi0,PCDstar,plotstyle{i},'LineWidth',3)
    leg{i} = sprintf('\\mu_0 = %d',mu0(i));
end
figure(1)
xlabel('\xi_0');ylabel(['Optimal trial length T^*_',subscript])
legend(leg,'Location','North')
hold off
UtilStdizeFigureAbs(1);
figure(2)
xlabel('\xi_0');ylabel(['Maximum expected net gain B^*_',subscript])
legend(leg,'Location','North')
hold off
UtilStdizeFigureAbs(2);
figure(3)
xlabel('\xi_0');ylabel(['Optimal number of pairwise allocations Q^*_',subscript])
legend(leg,'Location','North')
hold off
UtilStdizeFigureAbs(3);
figure(4)
xlabel('\xi_0');ylabel(['Probability of correct decision PCD_',subscript])
legend(leg,'Location','North')
hold off
UtilStdizeFigureAbs(4);

%%%%%%%%%%%%%%%%%%%%%
% Power curves
%%%%%%%%%%%%%%%%%%%%%
basic = MixInputModifier(basic,{'mu0',0});
xi0 = [2 5 20];
plotstyle = {'r:','b-','k--'};
leg = cell(size(xi0));
figure(5)
hold on
for i = 1:length(xi0)
    basic = MixInputModifier(basic,{'UnkVarianceShape',xi0(i)});
    [Tstar,~] = OneShotMaxExpectedNetGain(basic);
    [power,w] = OneShotCPCD(basic,Tstar);
    plot(w,power,plotstyle{i},'LineWidth',3)
    leg{i} = sprintf('\\xi_0 = %0.1f',xi0(i));
end
ylim([0 1])
legend(leg)
xlabel('True expected INMB, W')
ylabel('Power')
hold off
UtilStdizeFigureAbs(5);
basic = MixInputModifier(basic,{'UnkVarianceShape',3});

% Save the figures if dosave is true
if nargin < 3
    dosave = false;
end

if dosave
    UtilSaveFigFile(1, 'Profher', 'OptTrialLength_vs_varBeliefs',[subonline,horizon],'eps');
    UtilSaveFigFile(2, 'Profher', 'MaxENG_vs_varBeliefs',[subonline,horizon],'eps');
    UtilSaveFigFile(3, 'Profher', 'OptNumPairs_vs_varBeliefs',[subonline,horizon],'eps');
    UtilSaveFigFile(4, 'Profher', 'PCD_vs_varBeliefs',[subonline,horizon],'eps');
end

