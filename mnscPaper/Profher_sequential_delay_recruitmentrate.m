function Profher_sequential_delay_recruitmentrate(dosave)
% This function takes in the order of 1h
close all
[basic, advanced] = ProfherSetParameters_sequential();
[basic, advanced, rval, msgs] = DelayInputValidator( basic, advanced );
if ~rval, msgs, end
basicmix = ProfherSetParameters();


monthlydiscountrate = (1 + 0.035)^(1/12) - 1 ; % 
incrate = 7000/12; % incidence rate
nummonths = 15*12;   % number of months of validity of decision
basePopValue = incrate * nummonths;

delayinmonths = [1:11,12,13:18]; % delay in observing INMB in the study ( in months )  
patientspermonth = [3:7,94/12,9:12];  % maximum number of patients put into trial per month 
ENGseq = zeros(length(delayinmonths),length(patientspermonth));
ENGoneshot = ENGseq;
Qseq = ENGseq; % Expected number of samples
Qoneshot = ENGseq; % Number of samples for the oneshot design
fignum = 20;
tic
for i = 1:length(delayinmonths)
    for j = 1:length(patientspermonth)
        basic.tau = delayinmonths(i) * patientspermonth(j)/2 ;
        basic.theta = (1 + monthlydiscountrate)^(-2/patientspermonth(j)); %  per patient pair discount rate
        DTdiscpatient = (1 + monthlydiscountrate)^(1/patientspermonth(j)) - 1; % DT discount parameter per patient
        basic.P = (1 - (1 + DTdiscpatient)^(-patientspermonth(j)*basePopValue/incrate))/(1 - (1 + DTdiscpatient)^(-patientspermonth(j)/incrate)); % number of patients to benefit from the adoption decision
        % Start running the analysis
        [~, mat] = DelayCurvesRecur(basic, advanced);
        [mat] = DelayStageOne(basic, advanced, mat );
        toneshot = (0:basic.TMax);
        [mat]  = DelayOptimalBayesOneStage ( basic, advanced, toneshot, mat ); % this function that Steve wrote defines the optimal fixed sample size based on a comparison of EVSI and sampling costs
%         [ fignum, mat ] = DelaySimOverview( fignum, basic, advanced, mat ) ;
        [~,OneShotIndex] = min(abs(mat.muvec - basic.mu0));
        % B0hat instead B0vec? Yes, B0hat is value-to-go at Delta units,
        % and B0vec is value-to-go at 0 time units
        ENGseq(i,j) = mat.B0hat(OneShotIndex) - basicmix.ccap(patientspermonth(j));
%         ENGseq(i,j) = mat.B0vec(OneShotIndex) - basicmix.ccap(patientspermonth(j));
        ENGoneshot(i,j) = mat.OptOneShotReward(OneShotIndex) - basicmix.ccap(patientspermonth(j));
        Qseq(i,j) = mat.optENumSamps(OneShotIndex);
        Qoneshot(i,j) = mat.OptOneShotLength(OneShotIndex);
        disp([i,j])
    end
end
toc

figure(1)
plot(delayinmonths,ENGseq(:,6)./1e6,'k-',delayinmonths,ENGoneshot(:,6)./1e6,'b--','LineWidth',2)
hold on
xlabel('Delay (mo)');
ylabel('Expected net gain V (millions of £)')
legend('Sequential','One-shot')
UtilStdizeFigureAbs(1);

figure(2)
plot(delayinmonths,ENGseq(:,6)./ENGoneshot(:,6),'k-','LineWidth',2)
hold on
xlabel('Delay (mo)');
ylabel('Ratio of expected net gain (sequential to one-shot)')
UtilStdizeFigureAbs(2);

figure(3)
plot(patientspermonth,ENGseq(12,:)./1e6,'k-',patientspermonth,ENGoneshot(12,:)./1e6,'b--','LineWidth',2)
hold on
xlabel('Recruitment rate (mo^{-1})');
ylabel('Expected net gain V (millions of £)')
legend('Sequential','One-shot','Location','SouthEast')
UtilStdizeFigureAbs(3);

figure(4)
plot(patientspermonth,Qseq(12,:),'k-',patientspermonth,Qoneshot(12,:),'b--','LineWidth',2)
hold on
xlabel('Recruitment rate (mo^{-1})');
ylabel('Expected sample size Q')
legend('Sequential','One-shot','Location','SouthEast')
UtilStdizeFigureAbs(4);

% Save the figures and data if dosave=true
if nargin < 1
    dosave = false;
end
if dosave
    save('Profher_sequential_delay_recruitmentrate.mat','ENGseq','ENGoneshot','Qseq','Qoneshot','patientspermonth','delayinmonths');
    UtilSaveFigFile(1, 'Sequential', 'ProfherSequentialDelay',[],'eps');
    UtilSaveFigFile(2, 'Sequential', 'ProfherSequentialDelay_relative',[],'eps');
    UtilSaveFigFile(3, 'Sequential', 'ProfherSequentialRecRate',[],'eps');
    UtilSaveFigFile(4, 'Sequential', 'ProfherSequentialRecRate_Q',[],'eps');
end
