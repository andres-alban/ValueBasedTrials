function ProfherRecRate(horizon,online,dosave)
% This function evaluates the sensitivity to the recruitment rate. 
% It creates some of the figures presented in section 7 of the paper
% 
%  Revision: AA 29-11-2017

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




% This chunk plots a graph of Tstar and ENGstar vs r for a three
% different values of mu0 
zeta = basic.zeta;
r = 10:10:zeta;
mu0 = -1500:1500:1500;
plotstyle = {'r:','b-','k--'};
leg = cell(size(mu0));
for i = 1:length(mu0)
    basic = MixInputModifier(basic,{'mu0',mu0(i)});
    [Tstar,ENGstar,PCDstar,p] = SensitivityAnalysisOneParameter(basic,'r',r);
    figure(1)
    hold on
    plot(p,Tstar,plotstyle{i},'LineWidth',3)
    figure(2)
    hold on
    Qstar = Tstar .* p;
    plot(p,Qstar,plotstyle{i},'LineWidth',3)
    figure(3)
    hold on
    plot(p,ENGstar,plotstyle{i},'LineWidth',3)
    figure(4)
    hold on
    plot(p,PCDstar,plotstyle{i},'LineWidth',3)
    
    leg{i} = sprintf('\\mu_0 = %d',mu0(i));
end
figure(1)
xlabel('Recruitment rate');ylabel(['Optimal trial length T^*_',subscript])
legend(leg)
hold off
UtilStdizeFigureAbs(1);
figure(2)
xlabel('Recruitment rate');ylabel(['Number of pairwise allocations Q^*_',subscript])
legend(leg,'Location','NorthWest')
hold off
UtilStdizeFigureAbs(2);
figure(3)
xlabel('Recruitment Rate');ylabel(['Maximum expected net gain B^*_',subscript])
legend(leg)
hold off
UtilStdizeFigureAbs(3);
figure(4)
xlabel('Recruitment Rate');ylabel(['Probability of correct decision PCD_',subscript])
legend(leg,'Location','SouthEast')
hold off
UtilStdizeFigureAbs(4);


%%%%%%%%%%%%%%%%%%%%%
% Power curves
%%%%%%%%%%%%%%%%%%%%%
basic = MixInputModifier(basic,{'mu0',0});
r = [50,100,1000];
plotstyle = {'r:','b-','k--'};
leg = cell(size(r));
figure(5)
hold on
for i = 1:length(r)
    basic = MixInputModifier(basic,{'r',r(i)});
    [Tstar,~] = OneShotMaxExpectedNetGain(basic);
    [power,w] = OneShotCPCD(basic,Tstar);
    plot(w,power,plotstyle{i},'LineWidth',3)
    leg{i} = sprintf('r = %d',r(i));
end
ylim([0 1])
legend(leg)
xlabel('True expected INMB, W')
ylabel('Power')
hold off
UtilStdizeFigureAbs(5);
basic = MixInputModifier(basic,{'r',50});


        
% Save the figures if dosave=true
if nargin < 3
    dosave = false;
end
if dosave
    UtilSaveFigFile(1, 'Profher', 'OptTrialLength_vs_r',[subonline,horizon],'eps');
    UtilSaveFigFile(2, 'Profher', 'OptNumPairs_vs_r',[subonline,horizon],'eps');
    UtilSaveFigFile(3, 'Profher', 'MaxENG_vs_r',[subonline,horizon],'eps');
    UtilSaveFigFile(4, 'Profher', 'PCD_vs_r',[subonline,horizon],'eps');
    UtilSaveFigFile(5, 'Profher', 'CPCD_vs_r',[subonline,horizon],'eps');
end



