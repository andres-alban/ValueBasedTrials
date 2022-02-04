function ProfherPatentMktExcl(online,dosave)
% This function compares  market exclusivity to
% patent protection. It plots some of the figues in section 7 of the paper.
% 
% Revision: AA 26-11-2017

close all
if nargin < 1
    online = 1;
end
if online == 1
    subonline = 'on';
elseif online == 0
    subonline = 'off';
else
    error('online wrongly specified. Choose 1 for online or 0 for offline')
end
% Define the model with the appropriate parameters for Profher
basic = ProfherSetParameters();
basic = MixInputModifier(basic,{'online',online});

horizon = {'Patent','MktExcl'};
mu0 = -1500:1500:1500;
bargraphTstar = zeros(length(horizon),length(mu0));
bargraphENGstar = bargraphTstar;
bargraphPCDstar = bargraphTstar;
leg = cell(1,length(mu0));
for i = 1:length(mu0)
    basic = MixInputModifier(basic,{'mu0',mu0(i)});
    [Tstar,ENGstar,PCDstar,~] = SensitivityAnalysisOneParameter(basic,'horizon',horizon);
    Tstar(Tstar==0) = 0.05;
    bargraphTstar(:,i) = Tstar';
    bargraphENGstar(:,i) = ENGstar';
    bargraphPCDstar(:,i) = PCDstar';
    leg{i} = sprintf('\\mu_0 = %d',mu0(i));
end
horizon = {'Patent Protection','Market Exclusivity'};
figure(1)
bar(bargraphTstar);
set(gca,'XTickLabel',horizon)
ylabel(['Optimal trial length T^*_{',subonline,'}'])
legend(leg,'Location','NorthWest')
colormap('gray')
UtilStdizeFigureAbs(1);
figure(2)
bar(bargraphENGstar)
set(gca,'XTickLabel',horizon)
ylabel(['Maximum expected net gain B^*_{',subonline,'}'])
legend(leg,'Location','NorthWest')
colormap('gray')
UtilStdizeFigureAbs(2);
figure(3)
bargraphQstar = basic.r.*bargraphTstar;
bar(bargraphQstar);
set(gca,'XTickLabel',horizon)
ylabel(['Optimal number of pairwise allocations Q^*_{',subonline,'}'])
legend(leg,'Location','NorthWest')
colormap('gray')
UtilStdizeFigureAbs(3);
figure(4)
bar(bargraphPCDstar)
set(gca,'XTickLabel',horizon)
ylabel(['Probability of correct decision PCD_{',subonline,'}'])
legend(leg,'Location','South')
colormap('gray')
UtilStdizeFigureAbs(4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sensitivity of horizon as a function of recruitment rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
basic = MixInputModifier(basic,{'mu0',0,'horizon','Patent'});
zeta = basic.zeta;
r = 10:10:(zeta);
[TstarPatent,ENGstarPatent,~,~] = SensitivityAnalysisOneParameter(basic,'r',r);
basic = MixInputModifier(basic,{'horizon','MktExcl'});
[TstarMktExcl,ENGstarMaktExcl,~,p] = SensitivityAnalysisOneParameter(basic,'r',r);
figure(5)
plot(p,TstarPatent,'b-',p,TstarMktExcl,'r:','LineWidth',3)
xlabel('Recruitment rate');ylabel('Optimal trial length')
legend('Patent Protection','Market Exclusivity')
UtilStdizeFigureAbs(5);
figure(6)
plot(p,ENGstarPatent,'b-',p,ENGstarMaktExcl,'r:','LineWidth',3)
xlabel('Recruitment Rate');ylabel(['Maximum expected net gain B^*_{',subonline,'}'])
legend('Patent Protection','Market Exclusivity','Location','East')
UtilStdizeFigureAbs(6);

%%%%%%%%%%%%%%%%%%%%%
% Power curves
%%%%%%%%%%%%%%%%%%%%%
basic = MixInputModifier(basic,{'mu0',0,'horizon','Patent'});
horizon = {'Patent','MktExcl'};
plotstyle = {'r:','b-'};
figure(7)
hold on
for i = 1:length(horizon)
    basic = MixInputModifier(basic,{'horizon',horizon(i)});
    [Tstar,~] = OneShotMaxExpectedNetGain(basic);
    [power,w] = OneShotCPCD(basic,Tstar);
    plot(w,power,plotstyle{i},'LineWidth',3)
end
ylim([0 1])
legend(horizon)
xlabel('Given expected INMB, W')
ylabel('Power')
hold off
UtilStdizeFigureAbs(7);
basic = MixInputModifier(basic,{'horizon','Patent'});



% Plot figures if dosave is true
if nargin < 2
    dosave = false;
end

if dosave
    UtilSaveFigFile(1, 'Profher', 'OptTrialLength_vs_horizon',subonline,'eps');
    UtilSaveFigFile(2, 'Profher', 'MaxENG_vs_horizon',subonline,'eps');
    UtilSaveFigFile(3, 'Profher', 'OptNumPairs_vs_horizon',subonline,'eps');
    UtilSaveFigFile(4, 'Profher', 'PCD_vs_horizon',subonline,'eps');
    UtilSaveFigFile(5, 'Profher', 'OptTrialLength_horizon_r',subonline,'eps');
    UtilSaveFigFile(6, 'Profher', 'MaxENG_horizon_r',subonline,'eps');
    UtilSaveFigFile(7, 'Profher', 'CPCD_vs_horizon',subonline,'eps');
end