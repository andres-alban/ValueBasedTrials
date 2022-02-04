function ProfherSensitivityBeliefs(horizon,online,dosave)
% This function evaluates the sensitivity of the maximization problem to
% the beliefs mu0 and n0. It plots some of the figures in section 7 of the
% paper
% 
% Revision: AA 29-11-2017

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


% Plot 3D figures of the Tstar and ENGstar as a function of both mu0 and n0
mu0 = -1500:30:1500;
n0log = -1:0.1:2.5;
n0 = 10.^n0log;
[Tstar,ENGstar,PCDstar,p1,p2] = SensitivityAnalysisTwoParameters(basic,'mu0',mu0,'n0',n0);
figure(1)
step = 0.1;
range = 6*step;
v = round(max(max(Tstar))/step)*step;
v = (v-range):step:v;
[C,h] = contour(p1,p2,Tstar,v,'LineWidth',2);
set(gca,'YScale','log');
xlabel('\mu_0');ylabel('n_0');title(['Optimal trial length T^*_',subscript,' (yr)'])
clabel(C,h)
% view(142.5,30)
UtilStdizeFigureAbs(1);
figure(2)
ENGstar = ENGstar./1e6;
v = [1 2 3 4 7 10]*10;
[C,h] = contour(p1,p2,ENGstar,v,'LineWidth',2);
set(gca,'YScale','log');
xlabel('\mu_0');ylabel('n_0');title(['Maximum expected net gain B^*_',subscript,' (Mill. £)'])
clabel(C,h)
UtilStdizeFigureAbs(2);
figure(3)
step = 10;
range = 6*step;
Qstar = basic.r.*Tstar;
v = round(max(max(Qstar))/step)*step;
v = (v-range):step:v;
[C,h] = contour(p1,p2,Qstar,v,'LineWidth',2);
set(gca,'YScale','log');
xlabel('\mu_0');ylabel('n_0');title(['Optimal number of pairwise allocations Q^*_',subscript])
clabel(C,h)
UtilStdizeFigureAbs(3);
figure(6)
v = min(min(PCDstar));
v = round(v/0.1)*0.1;
v = v:0.05:1;
[C,h] = contour(p1,p2,PCDstar,v,'LineWidth',2);
set(gca,'YScale','log');
xlabel('\mu_0');ylabel('n_0');title(['Probability of correct decision PCD_',subscript])
clabel(C,h)
UtilStdizeFigureAbs(6);


% Plot Tstar and ENGstar as a function of mu0 only
[Tstar,ENGstar,~,p] = SensitivityAnalysisOneParameter(basic,'mu0',mu0);
figure(7)
plot(p,Tstar,'LineWidth',3)
xlabel('\mu_0');ylabel(['Optimal trial length T^*_',subscript])
UtilStdizeFigureAbs(7);
figure(8)
plot(p,ENGstar,'LineWidth',3)
xlabel('\mu_0');ylabel(['Maximum expected net gain B^*_',subscript])
UtilStdizeFigureAbs(8);

% Plot Tstar and ENGstar as a function of n0 only
[Tstar,ENGstar,~,p] = SensitivityAnalysisOneParameter(basic,'n0',n0);
figure(9)
semilogx(p,Tstar,'LineWidth',3)
xlabel('n_0');ylabel(['Optimal trial length T^*_',subscript])
UtilStdizeFigureAbs(9);
figure(10)
plot(p,ENGstar,'LineWidth',3)
xlabel('n_0');ylabel(['Maximum expected net gain B^*_',subscript])
UtilStdizeFigureAbs(10);

%%%%%%%%%%%%%%%%%%%%%
% Power curves
%%%%%%%%%%%%%%%%%%%%%
mu0 = [-1500,0,1500];
plotstyle = {'r:','b-','k--'};
leg = cell(size(mu0));
figure(11)
hold on
for i = 1:length(mu0)
    basic = MixInputModifier(basic,{'mu0',mu0(i)});
    [Tstar,~] = OneShotMaxExpectedNetGain(basic);
    [power] = OneShotCPCD(basic,Tstar);
    plot(w,power,plotstyle{i},'LineWidth',3)
    leg{i} = sprintf('\\mu_0 = %d',mu0(i));
end
ylim([0 1])
legend(leg)
xlabel('True expected INMB, W')
ylabel('Power')
hold off
UtilStdizeFigureAbs(11);
basic = MixInputModifier(basic,{'mu0',0});

%%%%%%%%%%%%%%%
n0 = [1/2,1,20];
plotstyle = {'r:','b-','k--'};
leg = cell(size(n0));
figure(12)
hold on
for i = 1:length(n0)
    basic = MixInputModifier(basic,{'n0',n0(i)});
    [Tstar,~] = OneShotMaxExpectedNetGain(basic);
    [power,w] = OneShotCPCD(basic,Tstar);
    plot(w,power,plotstyle{i},'LineWidth',3)
    leg{i} = sprintf('n_0 = %d',n0(i));
end
ylim([0 1])
legend(leg)
xlabel('True expected INMB, W')
ylabel('Power')
hold off
UtilStdizeFigureAbs(12);
basic = MixInputModifier(basic,{'n0',20});


% Save the figures if dosave is true
if nargin < 3
    dosave = false;
end

if dosave
    UtilSaveFigFile(1, 'Profher', 'OptTrialLengthvsBeliefs',[subonline,horizon],'eps');
    UtilSaveFigFile(2, 'Profher', 'MaxENGvsBeliefs',[subonline,horizon],'eps');
    UtilSaveFigFile(3, 'Profher', 'OptNumPairsvsBeliefs',[subonline,horizon],'eps');
    UtilSaveFigFile(6, 'Profher', 'PCDvsBeliefs',[subonline,horizon],'eps');
    UtilSaveFigFile(7, 'Profher', 'OptTrialLengthvsmu0',[subonline,horizon],'eps');
    UtilSaveFigFile(8, 'Profher', 'MaxENGvsmu0',[subonline,horizon],'eps');
    UtilSaveFigFile(9, 'Profher', 'OptTrialLengthvsn0',[subonline,horizon],'eps');
    UtilSaveFigFile(10, 'Profher', 'MaxENGvsn0',[subonline,horizon],'eps');
    UtilSaveFigFile(11, 'Profher', 'CPCDvsmu0',[subonline,horizon],'eps');
    UtilSaveFigFile(12, 'Profher', 'CPCDvsn0',[subonline,horizon],'eps');
end



