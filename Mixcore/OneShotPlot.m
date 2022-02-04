function OneShotPlot(basic,par,parVal,dosave,folder,file)
% This function plots the expected net gain as a function of the duration
% and the CPCS as a function of the true INMB for the trial defined by
% basic and an alternative that modifies parameter par to the value parVal.

basicTemp = basic;
linetype = {'b-','b-.'};
if nargin < 3
    par = [];
end
if isempty(par)
    K = 1;
else
    K = 2;
end

%% ENG plot
Tstar = zeros(1,K);
ENGstar = Tstar;
figure(1)
for i = 1:K
    if (i == 2)
        basicTemp = MixInputModifier(basicTemp,{par,parVal});
    end
    % Solve the maximization problem
    ENG = @(T) OneShotExpectedNetGain(basicTemp,T);
    [Tstar(i),ENGstar(i)] = OneShotMaxExpectedNetGain(basicTemp);
    % Calculate the sample size that would be prescribed by the standard power
    % calculations
    T = 0:(1/basicTemp.r):(basicTemp.Tmax);
    plot(T,ENG(T),linetype{i},'LineWidth',2);
    hold on
end
for i = 1:K
    plot(Tstar(i),ENGstar(i),'ro','LineWidth',2,'MarkerSize',6,'MarkerFaceColor','r')
end
N = FrequentistSampleSize(basicTemp);
Tfreq = 2*N/basicTemp.r;
hax = axis;
plot([Tfreq Tfreq],hax(3:4),'k--','LineWidth',2);

if K == 2
    legend('Base',[par,' = ',num2str(parVal)])
end
xlabel('Recruitment period length (mo)')
ylabel('Expected Net Gain (£)')
UtilStdizeFigureAbs(1);
hold off

%% CPCS plot

basicTemp = basic;
delta = basicTemp.relevantdifference;
w = linspace(-1.5.*delta,1.5*delta,delta);

figure(2)
CPCSbase = OneShotCPCD(w,basicTemp,Tstar(1));
plot(w,CPCSbase,linetype{1},'LineWidth',2)
hold on

if K == 2
basicTemp = MixInputModifier(basicTemp,{par,parVal});
CPCSalt = OneShotCPCD(w,basicTemp,Tstar(2));
plot(w,CPCSalt,linetype{2},'LineWidth',2)
end

plot([delta delta],[0 1],'k--',[-delta -delta],[0 1],'k--','LineWidth',1)
xlabel('w')
ylabel('CPCS')
if K==2
    legend('Base',[par,' = ',num2str(parVal)],'Location','SouthEast')
end
UtilStdizeFigureAbs(2);
hold off

%% Save the figure if dosave = true
if nargin < 4
    dosave = false;
end

if dosave
    if nargin < 5
        folder = 'Temp';
        file = folder;
    elseif nargin < 6
        file = folder;
    end
    UtilSaveFigFile(1, folder, file, ['_',par,'_',num2str(parVal)],'eps');
    UtilSaveFigFile(2, folder, file, ['_CPCS_',par,'_',num2str(parVal)],'eps');
end
