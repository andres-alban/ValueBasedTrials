function[ basic, advanced] = ProfherSetParameters_sequential()
% ProfherSetParameters_sequential: the input parameters for the ProFHER 
% sequential application section in the theory paper.
% Adapted from the JRSSB paper of Chick, Forster, Pertile
%(alpha order) and extended for the ProFHER theory paper.
%
% Created: 16 February 2021
% Last touched: 16 February 2021
% 
% 
%%%%%% Case study parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
incrate = 7000/12; % incidence rate
nummonths = 15*12;   % number of months of validity of decision

[basic, advanced] = DelayInputConstructor() ;
basic.pN = 0.39 ; % proportion using New at the start of the study
basic.sigma = 4400; % sampling sd for INMB (in £)
basic.t0 = 2 ; % effective sample size of the prior distribution ( must be greater than 1 )
delayinmonths = 12; % delay in observing INMB in the study ( in months )   
basic.numpairs = 125 ; % number of pairwise allocations in the study ( NB this is NOT TMax )
accrualperiod =  32 ; % actual duration of the recruitment period of trial ( in months )             
basic.ICostN = 0 ; % cost of switching all patients to New ( in £ )
basic.ICostS = 0 ; % cost of switching all patients to Standard ( in £ )
cResearch = 4080; %3200 ; % research cost per pairwise allocation ( in £ )
monthlydiscountrate = (1 + 0.035)^(1/12) - 1 ; % 
patientspermonth = 94/12 ;  % maximum number of patients put into trial per month
basic.tau = delayinmonths * patientspermonth/2 ;
basic.TMax = 200 ; % assumption 
%basic.mumax = 10 * 50 * 250000 / (basic.t0 + basic.TMax) ; % 500000 was the max s_Y in Martin-Paolo code for the same application
basic.mumax = 2e4 ; 
basic.mumin = -2e4 ; 
basic.mu0 = 0;
basic.z  = 0 ; % critical value for classical statistical inference test ( e.g. 1.96 )
basic.theta = (1 + monthlydiscountrate)^(-2/patientspermonth); %  per patient pair discount rate

basePopValue = incrate * nummonths;
if monthlydiscountrate == 0
    basic.P = basePopValue; % number of patients to benefit from the adoption decision
else
    DTdiscpatient = (1 + monthlydiscountrate)^(1/patientspermonth) - 1; % DT discount parameter per patient
    basic.P = (1 - (1 + DTdiscpatient)^(-patientspermonth*basePopValue/incrate))/(1 - (1 + DTdiscpatient)^(-patientspermonth/incrate)); % number of patients to benefit from the adoption decision
end

basic.factor = 1.0 ; % adjustment factor in denominator ( refer to Theory paper )
basic.c = cResearch ; % adjusted research cost per pairwise allocation
basic.PPatients = basic.P ; % adjusted patients to benefit
basic.ICost = basic.ICostN; % adjusted switching cost
basic.addon = 0;        % adjusted add-on at end of trial


%%%%%%%%%
basic.online = false;  % set to true if per patient reward is gained ( online learning )
advanced.fixedP = true;  % set to true if there are exactly PPatients to follow once a treatment or alternative is selected
        % set to false if the patients which are not included in trial due
        % to early stopping are also included in the treatment arm, that
        % is, if fixedP is false, then there are P + (TMAX - T) patients
        % which are included, where T is the stopping time (number of
        % patients) included in the trial (and T \in 0, 1, ..., TMAX).
advanced.nochangeallowed = false;    % set to false (the better default) meaning that it is ok
        % to change what you think is best after the added tau data points
        % come in. set this to true if you want to ignore the added tau
        % pending data points to make a decision.

%%%%%% Numerical computation parameters %%%%%%

%advanced.MinGridPerStdev = (basic.TMax-basic.t0)/s_step_t;    % if this is negative, code tries to use its negative as the delta-t
%advanced.MinGridPerStdev = 200; 
advanced.MinGridPerStdev = 20 ;    % minimum number of grid points per stdev of output, or negative to use dt=1
% For plotting
advanced.smoothed = true ;   % true or false if curve smoothing is desired.
advanced.DOPLOT = false ; % put true or false depending if you do or do not want dynamic plots
advanced.NumGridsForContours = 100; % number of time values for generating contour plots (also num plots in animation to use during recursion over t)
advanced.bigfontsize=16 ;
advanced.smallfontsize=14 ;
advanced.fontname='Helvetica';
advanced.titlestring='ProFHER';


advanced.dirstring='Sequential';

advanced.filestring='ProFHER';
advanced.PLOTSIMS = false ;

advanced.saveplot = true;           % set to true to save plots in files, false to not save files automatically

end

