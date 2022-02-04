function ProfherSensitivitySwitchingCost(horizon,online,dosave)
% This function evaluates the sensitivity of the maximization problem to
% the beliefs mu0 and n0. It plots some of the figures in section 7 of the
% paper
% 
% Revision: AA 6-11-2017

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
IN = 0:1e6:1e8;
IS = 0:1e6:1e8;
[Tstar,ENGstar,PCDstar,p1,p2] = SensitivityAnalysisTwoParameters(basic,'IN',IN,'IS',IS);
figure(1)
[C,h] = contour(p1,p2,Tstar,5);
xlabel('I_N');ylabel('I_S');title('Optimal trial length')
clabel(C,h)
UtilStdizeFigureAbs(1);
figure(2)
[C,h] = contour(p1,p2,ENGstar,5);
xlabel('I_N');ylabel('I_S');title('Maximum expected net gain')
clabel(C,h)
UtilStdizeFigureAbs(2);
figure(3)
Qstar = basic.r.*Tstar;
[C,h] = contour(p1,p2,Qstar,5);
xlabel('I_N');ylabel('I_S');title('Optimal number of pairwise allocations')
clabel(C,h)
UtilStdizeFigureAbs(3);
figure(4)
[C,h] = contour(p1,p2,PCDstar,5);
xlabel('I_N');ylabel('I_S');title('Probability of correct decision')
clabel(C,h)
UtilStdizeFigureAbs(4);

%%%%%%%%%%%%%%%%%%%%%
% Power curves
%%%%%%%%%%%%%%%%%%%%%
basic = MixInputModifier(basic,{'mu0',0});
IN = [0 0 1e7 1e7];
IS = [0 1e7 0 1e7];
plotstyle = {'r:','b-','k--','m-.'};
leg = cell(size(IN));
figure(5)
hold on
for i = 1:length(IN)
    basic = MixInputModifier(basic,{'IN',IN(i),'IS',IS(i)});
    [Tstar,~] = OneShotMaxExpectedNetGain(basic);
    [power,w] = OneShotCPCD(basic,Tstar);
    plot(w,power,plotstyle{i},'LineWidth',3)
    leg{i} = sprintf('IN = %d, IS = %d',IN(i),IS(i));
end
ylim([0 1])
legend(leg)
xlabel('True expected INMB, W')
ylabel('Power')
hold off
UtilStdizeFigureAbs(5);


% Save the figures if dosave is true
if nargin <3
    dosave = false;
end

if dosave
    UtilSaveFigFile(1, 'Profher', 'OptTrialLengthvsI',[subonline,horizon],'eps');
    UtilSaveFigFile(2, 'Profher', 'MaxENGvsI',[subonline,horizon],'eps');
    UtilSaveFigFile(3, 'Profher', 'OptNumPairsvsI',[subonline,horizon],'eps');
    UtilSaveFigFile(4, 'Profher', 'PCDvsI',[subonline,horizon],'eps');
    UtilSaveFigFile(5, 'Profher', 'CPCDvsI',[subonline,horizon],'eps');
end
